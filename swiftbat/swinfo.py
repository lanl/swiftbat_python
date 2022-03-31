#! /usr/bin/env python


from __future__ import print_function, division, absolute_import

import functools
from pathlib import Path

"""
swinfo
More utilities for dealing with Swift Data
David Palmer   palmer@lanl.gov


"""

"""
Copyright (c) 2018, Triad National Security, LLC. All rights reserved.
 
This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL),
which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration.

All rights in the program are reserved by Triad National
Security, LLC, and the U.S. Department of Energy/National
Nuclear Security Administration. The Government is granted for
itself and others acting on its behalf a nonexclusive, paid-up,
irrevocable worldwide license in this material to reproduce,
prepare derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit others to
do so.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1.       Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2.       Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.       Neither the name of Triad National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import os
import sys

# add the path of this module to the searchpath to let helpers in
execdir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
if not execdir in map(os.path.abspath, sys.path):
    sys.path.append(execdir)

import re
import datetime
import shutil

import glob
import gzip
import traceback
import getopt
import time
import math
from swiftbat import swutil
from swiftbat.clockinfo import utcf
import astropy.units as u
from astropy.io import fits
from astropy import coordinates
import numpy as np
from io import StringIO

# Running into certificate problems from 2018-10-23 for
# https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error

import ssl

unsafe_context = ssl._create_unverified_context()

from urllib.request import urlopen, quote
from swiftbat import generaldir
import bs4 as BeautifulSoup
import ephem

import sqlite3  # Good sqlite3 tutorial at http://www.w3schools.com/sql/default.asp
split_translate = str.maketrans("][, \t;", "      ")



# Pointings database (as opposed to an observations database) has two tables, 'pointings' and 'days'
# pdbfile = os.path.join(swiftcache.theCache.params['LOCAL'][0],"pointings.db")
pdbfile = '/Volumes/Data/Swift/swift/pointings.db'
pdbfields_text = '''time_begin text unique primary key, time_end text, seqnum text, target_name text, ra_pnt real, dec_pnt real, roll_pnt real, asflown boolean'''
pdbfields = [s.split() for s in pdbfields_text.split(',')]
pdbfieldnames = [s[0] for s in pdbfields]

pdbdayfields_text = '''date text unique primary key, asflown_complete boolean, asflow_partial boolean, preplanned boolean'''
pdbdayfields = [s.split() for s in pdbdayfields_text.split(',')]
pdbdayfieldnames = [s[0] for s in pdbdayfields]


def adapt_boolean(bol):
    if bol:
        return "True"
    else:
        return "False"


def convert_boolean(bolStr):
    if str(bolStr) == "True":
        return bool(True)
    elif str(bolStr) == "False":
        return bool(False)
    else:
        raise ValueError("Unknown value of bool attribute '%s'" % bolStr)


sqlite3.register_adapter(bool, adapt_boolean)
sqlite3.register_converter("boolean", convert_boolean)

basecatalog = os.path.join(execdir, "catalog")
fitscatalog = os.path.join(execdir, "recent_bcttb.fits.gz")
# FIXME should be handled by dotswift
catalogdir = "/opt/data/Swift/analysis/sourcelist"
newcatalogfile = os.path.join(catalogdir, "newcatalog")
cataloglist = [os.path.join(catalogdir, "trumpcatalog"),
               os.path.join(catalogdir, "grbcatalog"),
               os.path.join(catalogdir, "catalog"),
               newcatalogfile]

earthradius_m = 6378100.0  # meters
r2d = 180 / 3.1415926  # this should be somewhere
atm_thickness = 100e3  # How deep in the atmosphere you should report attenuation

# define machineReadable=False


asflownURLpattern = "http://www.swift.psu.edu/operations/obsSchedule.php?d=%Y-%m-%d&a=1"
preplannedURLpattern = "http://www.swift.psu.edu/operations/obsSchedule.php?d=%Y-%m-%d&a=0"

ydhms = "%Y-%j-%H:%M:%S"

# TLEpattern = ["ftp://heasarc.gsfc.nasa.gov/swift/data/obs/%Y_%m/",".*","auxil","SWIFT_TLE_ARCHIVE.*.gz"]
TLEpattern = ["http://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/%Y_%m/", ".*", "auxil", "SWIFT_TLE_ARCHIVE.*.gz"]

tlefile = "/tmp/latest_swift_tle.gz"
tlebackup = os.path.expanduser("~/.swift/recent_swift_tle.gz")

radecnameRE = re.compile(r'''^(?P<rastring>[0-9.]+)_+(?P<decstring>([-+]*[0-9.]+))$''')

ipntimeRE = re.compile(r"""'(?P<mission>[^']+)'\s*'(?P<d>[ 0-9]+)/(?P<m>[ 0-9]+)/(?P<y>[ 0-9]+)'\s+(?P<s>[0-9.]+)""")

_maxobslength = datetime.timedelta(minutes=90)

verbose = False  # Verbose controls diagnositc output
terse = False  # Terse controls the format of ordinary output


# Old fashioned
def simbadnames(query):
    """Given a source name, or other SIMBAD query, generates a list of identifier matches
       See http://simbad.u-strasbg.fr/simbad/sim-fscript
    """
    u = urlopen(
        """http://simbad.u-strasbg.fr/simbad/sim-script?submit=submit+script&script=format+object+%%22+%%25IDLIST%%22%%0D%%0Aquery++%s""" % (
            quote(query),), None, 60)
    names = []
    while u:
        l = u.readline().decode()
        if re.match("""^::error::*""", l):
            raise ValueError(query)
        if re.match("""^::data::*""", l):
            break
    for l in u.readlines():
        s = l.decode().strip().split()
        if len(s) > 0:
            if s[0] == 'NAME':
                s = s[1:]
            names.append(" ".join(s))
    return names


# Astroquery
def simbadlocation(objectname):
    """
    Location according to 
    :param objectname: 
    :return: (ra_deg, dec_deg)
    """
    from astroquery.simbad import Simbad
    try:
        table = Simbad.query_object(objectname)
        if len(table) != 1:
            raise RuntimeError(f"No unique match for {objectname}")
        co = coordinates.SkyCoord(table['RA'][0], table['DEC'][0],
                                  unit=(u.hour, u.deg), frame='fk5')
        return (co.ra.degree, co.dec.degree)
    except Exception as e:
        raise RuntimeError(f"{e}")


class orbit:
    def __init__(self):
        self.getTLE()

    def getTLE(self):
        # TLE in this case is Two Line Element, not 3
        global verbose
        if verbose:
            print("Updating TLE")
        self.updateTLE()
        if verbose:
            print("Reading TLE from %s" % (tlefile))
        allTLE = gzip.open(tlefile, "rt" if sys.version_info.major == 3 else "r").readlines()
        # Some TLE files seem to have additional blank lines
        # Actually, what they have is \r\n in some cases and \n in others
        # Nevertheless, put this in in case they change the format
        allTLE = [l.strip() for l in allTLE if len(l.strip()) == 69]
        nTLE = len(allTLE) // 2
        self._tleByTime = [(datetime.datetime(2000 + int(allTLE[2 * i][18:20]), 1, 1, 0, 0, 0)
                            + datetime.timedelta(float(allTLE[2 * i][20:32]) - 1),
                            ("Swift", allTLE[2 * i], allTLE[2 * i + 1]))
                           for i in range(nTLE)]
        # for (ttle,tle) in self._tleByTime :
        #    print ("%s" %ttle)
        self.pickTLE(datetime.datetime.utcnow())
        self._earthcenter = ephem.Observer()
        self._earthcenter.lat = 0.0
        self._earthcenter.long = 0.0
        self._earthcenter.elev = -6378137.0  # equatorial radius of WGS84
        self._earthcenter.date = ephem.now()

    def pickTLE(self, t):
        if self._tleByTime[-1][0] < t:
            self._tle = self._tleByTime[-1][1]
            self._tletimes = [self._tleByTime[-1][0], datetime.datetime(datetime.MAXYEAR, 1, 1)]
        else:
            for w in range(len(self._tleByTime) - 1, -1, -1):
                if self._tleByTime[w][0] < t or w == 0:
                    self._tle = self._tleByTime[w][1]
                    self._tletimes = [self._tleByTime[w][0], self._tleByTime[w + 1][0]]
                    break

    def getSatellite(self, t):
        global verbose
        # print(t, self._tletimes)
        if t < self._tletimes[0] or self._tletimes[1] < t:
            if verbose:
                print("Picking TLE for time %s" % (t,))
            self.pickTLE(t)
        if verbose:
            print("Picked TLE for %s < %s < %s" % (self._tletimes[0], t, self._tletimes[1]))
        sat = ephem.readtle(self._tle[0], self._tle[1], self._tle[2])
        self._earthcenter.date = t
        sat.compute(self._earthcenter)
        return sat

    def usetledb(self, catnum=28485):
        """Use spacetrack tles from database"""
        import tledb
        tles = tledb.getTLEs(catnums=[catnum], all_in_time=True)
        self._tleByTime = [(t.epoch, t.threelines(split=True)) for t in tles]
        return self.getSatellite(self._earthcenter.date.datetime())

    def satname(self):
        return self._tleByTime[0][1][0]

    def updateTLE(self):
        global verbose
        try:
            # time.clock() is not what was wanted, since that doesn't give you the clock time
            checksecs = time.mktime((datetime.datetime.now() + datetime.timedelta(days=-1)).timetuple())
            if os.stat(tlefile).st_mtime > checksecs:
                return  # The TLE file exists and is less than a day old
        except:

            pass
        tlematch = re.compile(TLEpattern[-1])
        for url_month in (datetime.datetime.utcnow().strftime(TLEpattern[0]),
                          (datetime.datetime.utcnow() - datetime.timedelta(30)).strftime(TLEpattern[0])):
            try:
                if verbose:
                    print("Trying to get TLE from %s" % (url_month,))
                httpdir = generaldir.httpDir(url_month)
                obsmonth = httpdir.dirs()
                obsmonth.reverse()
                for obs in obsmonth:
                    auxfiles = httpdir.files(obs + "/auxil")
                    for f in auxfiles:
                        if tlematch.match(f):
                            if verbose:
                                print("Copying TLEs from " + obs + "/auxil/" + f)
                            httpdir.copyToFile(obs + "/auxil/" + f, tlefile)
                            try:
                                shutil.copyfile(tlefile, tlebackup)
                            except:
                                pass
                            return
            except Exception as e:
                pass
        try:
            if verbose:
                print("Using old verion of Swift TLE from " + tlebackup)
            shutil.copyfile(tlebackup, tlefile)
            shutil.copystat(tlebackup, tlefile)
            return
        except:
            raise LookupError("Could not find a recent TLE file")


# Source, initialized from a data string from the catalog files


def batExposure(theta, phi):
    """
    Given theta,phi in radians, returns (open_coded_area_in_cm^2, cosfactor)

    theta = distance from center of FOV (boresight) in radians
    phi = angle around boresight in radians.
    This is moderate accuracy, but does not take into account, e.g. dead detectors.
    """
    if math.cos(theta) < 0:
        return (0.0, math.cos(theta))
    # BAT dimensions
    detL = 286 * 4.2e-3  ## Amynote: each det element is has length 4.2e-3 m
    detW = 173 * 4.2e-3  ## Amynote: each det element is has length 4.2e-3 m
    maskL = 487 * 5.0e-3  # Using the 5 mm mask cells, true size is 95.9" x 47.8"
    maskW = 243 * 5.0e-3
    efl = 1.00
    # Calculate dx, but reverse it if necessary to put projected detector to left of mask center
    dx = (maskL - detL) / 2 - efl * math.tan(theta) * abs(math.cos(phi))
    dy = (maskW - detW) / 2 + efl * math.tan(theta) * math.sin(phi)
    # Boundaries of the detector as clipped to rectangle of mask
    x1 = max(0, -dx)
    x2 = min(detL, maskL - dx)
    deltaX = x2 - x1
    y1 = max(0, -dy)
    y2 = min(detW, maskW - dy)
    deltaY = y2 - y1
    # Now adjust for the cut corners, which first attacks the upper left detector corner
    # as described by the (if-necessary-reversed) coord. system
    xint = (y1 + dy) - maskW / 2 - (dx + x1)  # line of corner clip extends through (x1 + xint, y1)
    if deltaX < 0 or deltaY < 0:
        area = 0
    elif xint <= - deltaY:  # no clipping
        area = deltaX * deltaY
    elif xint <= 0:  # upper left corner clipped along left edge
        if deltaY <= deltaX - xint:  # clip from left edge, to top edge
            area = deltaX * deltaY - ((deltaY + xint) ** 2) / 2
        else:  # clip from left edge to right edge
            area = deltaX * -xint + (deltaX ** 2) / 2
    elif xint <= deltaX:  # clipped at bottom edge
        if xint <= deltaX - deltaY:  # clipped bottom and top
            area = (deltaX - xint) * deltaY - (deltaY ** 2) / 2
        else:  # clipped bottom and right
            area = ((deltaX - xint) ** 2) / 2
    else:
        area = 0
    # if you want to see what the corners do: area = max(0,deltaX * deltaY) - area
    # multiply by 1e4 for cm^2, 1/2 for open area
    return (area * 1e4 / 2, math.cos(theta))


def detid2xy(detids):
    """
    Convert detector ids to x,y positions

    This is tricky.  You can understand it, but it isn't worth it.

    :param detids: detector ids (ints 0...32767)
    :type detids: scalar int or convertible to a numpy array
    :return: x,y
    :rtype: x and y are each numpy arrays or scalars the same size as the original detids
    """
    scalar = np.isscalar(detids)
    detids = np.asarray(detids, dtype=np.int16)

    #  'mod' is a 128-detector half-DM
    block_and_mod, det_in_mod = np.divmod(detids, np.int16(128))
    block, mod_in_block = np.divmod(block_and_mod, np.int16(16))

    # x in {0..285}, y in {0..172}
    # origin lower left is block 8, hdm 9, det 0
    # 16 modules across (each 16 detectors wide) with 15 gaps of 2 -> 286 x values
    # 16 modules high (each 8 detectors high) with gaps of 3 -> 173 y values

    # Module frame imod, jmod
    # 127 119 111 103  95  87  79  71  56  48  40  32  24  16  8   0
    # 126                          70  57                          1
    # 125                          69  58                          2   ^
    # 124                          68  59                          3   |
    # 123                          67  60                          4   jmod
    # 122                          66  61                          5
    # 121                          65  62                          6    imod  ->
    # 120 112 104  96  88  80  72  64  63  55  47  39  31  23  15  7
    bit6, bit5_0 = np.divmod(det_in_mod, np.int16(64))  # bit7 is left or right half
    bit6_3, bit2_0 = np.divmod(det_in_mod, np.int16(8))  # bits 6-3 determine the imod, bits 2-0 jmod
    imod = np.int16(15) - bit6_3
    jmod = np.where(bit6, bit2_0, np.int16(7) - bit2_0)

    # for block 0-7, half-DMs are arranged
    # 1    9
    # 0    8
    # 3    11     ^
    # 2    10     |
    # 5    13     |
    # 4    12   jblock
    # 7    15
    # 6    14          iblock ---->

    bit10, bit9_7 = np.divmod(mod_in_block, np.int16(8))
    iblock = imod + np.int16(18) * bit10  # 16 detectors + 2-space gap for second column
    jmodrow = (np.array([6, 7, 4, 5, 2, 3, 0, 1]) * 11).astype(
        np.int16)  # 8 detectors and 3-space gap for each row of modules
    jblock = jmod + jmodrow[bit9_7]

    # blocks are arranged:
    #       0  1  2  3  4  5  6  7
    #       8  9  10 11 12 13 14 15
    # and blocks 8-15 are rotated 180 degrees
    bit15, bit14_11 = np.divmod(block, np.int16(8))
    #                      8-15        0-7
    y = np.where(bit15, np.int16(85) - jblock, np.int16(88) + jblock)
    x = np.where(bit15, np.int16(34) - iblock, iblock) + np.int16(36) * bit14_11
    if scalar:
        x = np.int16(x)
        y = np.int16(y)
    return x, y

@functools.lru_cache(maxsize=0)
def _xy2detidmap():
    """
    Produce a detector map filled with detector IDs

    -1 for unpopulated detector locations

    :return: detids[y, x]
    :rtype: uint16, shape = (173, 286)
    """
    detids = np.arange(0,2**15)
    x,y = detid2xy(detids)
    result = np.full((y.max()+1, x.max()+1), np.int16(-1))
    result[y,x] = detids
    return result

def xy2detid(x, y):
    dmap = _xy2detidmap()
    return dmap[y, x]


def loadsourcecat():
    """Read in source catalog
    
    FIXME, should be more automatic
    """
    global sourcecat
    global verbose
    try:
        id(sourcecat)
    except NameError:
        if os.path.exists(cataloglist[3]):
            if verbose:
                print("Loading catalogs:\n %s" % "\n ".join(cataloglist[0:3]))
            sourcecat = sourcelist(cataloglist[0:3] + [fitscatalog], verbose=verbose)  # Exclude newcatalog
        else:
            if verbose:
                print("Using old catalog %s" % (basecatalog,))
            sourcecat = sourcelist((basecatalog, fitscatalog), verbose=verbose)
        if verbose:
            print("Loaded")
        # print(" ".join(sourcecat.allsources.keys()))


class source:
    def __init__(self, initstring=None, **kwargs):
        """
        A source location with the ability to calculate BAT-relative angles and exposure

        :param initstring: String from source table, or an ephem function name (such as 'Sun')
        :param kwargs: {ra:ra_deg, dec:dec_deg, <name:'a name'>, <catnum:catnum>}
        """
        if initstring in ephem.__dict__:
            self.needs_computing = True
            self.computable = ephem.__dict__[initstring]()
        elif 'ra' in kwargs and 'dec' in kwargs:
            if initstring:
                raise RuntimeError("Give ra=,dec= or an initstring but not both")
            self.eq = ephem.Equatorial(kwargs['ra'] * ephem.degree,
                                       kwargs['dec'] * ephem.degree)
            self.name = kwargs.get('name', 'unnamed')
            self.catnum = kwargs.get('catnum', 0)
        else:
            self.needs_computing = False
            s = initstring.split("|")
            self.name = s[3].strip().replace(" ", "_")
            self.catnum = int(s[5])
            self.eq = ephem.Equatorial(ephem.degrees(float(s[9]) * ephem.degree),
                                       ephem.degrees(float(s[10]) * ephem.degree))
        # All input is in degrees, output is in ephem angles

    @classmethod
    def source(cls, ra, dec, name="anonymous", catnum=-1):
        return cls("|||{}||{}||||{}|{}".format(name, catnum, ra, dec))

    def exposure(self, ra, dec, roll):
        # returns (projected_area, cos(theta))
        (thetangle, phi) = self.thetangle_phi(ra, dec, roll)
        if thetangle.norm >= ephem.halfpi:
            return 0.0, 0.0
        return batExposure(thetangle, phi)

    def distance(self, ra, dec):
        return ephem.separation((self.eq.ra, self.eq.dec),
                                (ephem.degrees(ra * ephem.degree), ephem.degrees(dec * ephem.degree)))

    def posang_from(self, ra, dec):
        """ Position angle East of North from the given location to self """
        # Stolen from idl posang.pro, with self as point 2
        decrad = ephem.degrees(dec * ephem.degree)
        radiff = self.eq.ra - ephem.degrees(ra * ephem.degree)
        ang = ephem.degrees(math.atan2(math.sin(radiff),
                                       math.cos(decrad) * math.tan(self.eq.dec) - math.sin(decrad) * math.cos(radiff)))
        return ang

    def thetangle_phi(self, ra, dec, roll):
        """ Source position in instrument FOV given instrument pointing direction
            returns (thetangle,phi) where thetangle is the angular distance from
            the boresight and phi is the angle from the phi=0 axis of BAT.
            
            theta = tan(thetangle) gives the theta we use, which is the projected distance to a flat plane,
            but it is not useful for thetangle > 90 degrees
        """
        thetangle = self.distance(ra, dec)
        # roll is 'roll left' (CCW as seen from behind looking out along boresight),
        # posang is CCW, phi is CCW from +Y so (posang - roll - 90deg) gives phi
        phi = self.posang_from(ra, dec) - ephem.degrees(roll * ephem.degree) - ephem.halfpi
        # znorm, in range (-pi, pi) is how it is reported on-board
        phi = (phi + ephem.pi) % ephem.twopi - ephem.pi
        # print("posang = %f E of N, roll = %f" % (self.posang_from(ra,dec),0+ephem.degrees( roll * ephem.degree)))
        return (thetangle, phi)

    def compute(self, t):
        if self.needs_computing:
            obs = ephem.Observer()
            obs.elev = -6378137.0  # equatorial radius of WGS84
            obs.date = t
            # print(obs, self.computable)
            self.computable.compute(obs)
            self.name = self.computable.name
            self.catname = 0
            # print(self.computable)
            self.eq = ephem.Equatorial(self.computable)

    def __str__(self):
        return "|||%s||%d||||%f|%f|" % (self.name, self.catnum, self.eq.ra / ephem.degree, self.eq.dec / ephem.degree)


# This is a truncated version of catalog.  It should become the parent of catalog someday (FIXME 2008-09-17)
class sourcelist:
    def __init__(self, filelist=None,
                 verbose=False):
        if filelist is None:
            filelist = [os.path.join(catalogdir, "catalog"), os.path.join(catalogdir, "grbcatalog")]
        self.allsources = {}
        for file in filelist:
            # print(file)
            try:
                self.addFromFile(file, verbose=verbose)
            except:
                traceback.print_exc()
                pass

    def addFromFile(self, file, verbose=False):
        if verbose:
            print(file)
        if ".fits" in file:
            for row in fits.getdata(file):
                aSource = source.source(ra=row.field('RA_OBJ'),
                                        dec=row.field('DEC_OBJ'),
                                        name=row.field('NAME'))
                if verbose:
                    print(aSource)
                # print(line,aSource.name)
                self.allsources[self.canonName(aSource.name)] = aSource
        else:
            for line in open(file).readlines():
                try:
                    aSource = source(line)
                    if verbose:
                        print(aSource)
                    # print(line,aSource.name)
                    self.allsources[self.canonName(aSource.name)] = aSource
                except:
                    pass

    def byName(self, name):
        try:
            return self.allsources[self.canonName(name)]
        except:
            raise NameError

    def byCatnum(self, c):
        for s in self.allsources:
            if s.catnum == c:
                return s
        return None

    def canonName(self, name):
        return name.lower().replace(" ", "_")[0:15].replace("_", "")


def strip_nbsp(s):
    """ Replace '&nbsp' with spaces and then strip (beginning and end) spaces.
        An interior &nbsp; will remain an interior space"""
    return re.sub("&nbsp;", " ", s).strip()


class pointingEntry:
    def __init__(self, tstart, tslewend, tend, ra, dec, roll, obs, segment, sourcename, planned):
        self._tstart = tstart
        self._tslewend = tslewend
        if tend == None:  # How to handle the stubs  FIXME doesn't handle the case of reading day n then n-1
            self._tend = tstart + _maxobslength
        else:
            self._tend = tend
        self._ra = ra
        self._dec = dec
        self._roll = roll
        self._obs = obs
        self._segment = segment
        self._sourcename = sourcename
        self._planned = planned

    def fromTableRow(row, preplanned):
        """ row is a table row parsed by BeautifulSoup in the format listed below
            Returns a pointingEntry if it can, otherwise raises a ValueError"""
        try:
            entries = row.findAll('td')

            tstart = datetime.datetime.strptime(strip_nbsp(soupCatStrings(entries[0], True)), "%Y-%m-%d %H:%M:%S")
            tend = datetime.datetime.strptime(strip_nbsp(soupCatStrings(entries[1], True)), "%Y-%m-%d %H:%M:%S")
            obs = int(strip_nbsp(soupCatStrings(entries[2], True)))
            segment = int(strip_nbsp(soupCatStrings(entries[3], True)))
            sourcename = strip_nbsp(soupCatStrings(entries[4], True))
            ra = float(strip_nbsp(soupCatStrings(entries[5], True)))
            dec = float(strip_nbsp(soupCatStrings(entries[6], True)))
            roll = float(strip_nbsp(soupCatStrings(entries[7], True)))
            # Remaining fields not relevant
            pe = pointingEntry(tstart, tstart, tend, ra, dec, roll, obs, segment, sourcename, preplanned)
            if verbose:
                print(pe)
            return pe
        except:
            if verbose:
                print("Failed to parse row", row)
            raise ValueError("pointingEntry: Failed to parse row: %s", row)

    fromTableRow = staticmethod(fromTableRow)

    def fContains(self, t):
        return self._tstart <= t and t < self._tend

    def compareTime(self, t):
        if t < self._tstart:
            return -2  # Before the observation
        elif t < self._tslewend:
            return -1  # Before stable
        elif t < self._tend:
            return 0  # During stable observation
        else:
            return 1  # after observation

    def fOverlaps(self, trange, tend=None):
        if tend != None and len(trange) == 1:
            return self.fOverlaps((trange, tend))  # Convert from fOverlaps(x,y) to fOverlaps( (x,y) )
        else:
            return self._tstart <= trange[1] and self._tend >= trange[0]

    def __str__(self):
        if self._planned:
            if self._planned == True:
                plannedstring = "(Preplanned)"
            else:
                plannedstring = "(planned %s)" % self._planned
        else:
            plannedstring = ""
        return "%08d%03d  %s - %s - %s [%7.3f, %7.3f, %7.3f] %s %s" % (
            self._obs, self._segment,
            self._tstart.strftime("%Y-%j-%H:%M:%S"),
            self._tslewend.strftime("%H:%M:%S"),
            self._tend.strftime("%H:%M:%S"),
            self._ra, self._dec, self._roll,
            self._sourcename,
            plannedstring)


# The pre-planned is like
# 2007-086-00:00:00 | PPT | Begin | J2043+1255 | 36281 | 2 | 33590713 | 310.817427802637 | 12.9223455711058 | 85.6387273417218 | 0 | 0 | 12525 | 50 | 1471 | Paolo Giommi and Roger Romani Fill in Target | targ_20070850000_20070870000_01_DB.tcl
# 2007-086-00:00:00 | mnv | Begin | J2043+1255
# 2007-086-00:02:00 | mnv | End | J2043+1255
# 2007-086-00:27:00 | PPT | End | J2043+1255 | 36281 | 2 | 33590713 | 310.817427802637 | 12.9223455711058 | 85.6387273417218 | 0 | 0 | 12525 | 50 | 1471 | Paolo Giommi and Roger Romani Fill in Target | targ_20070850000_20070870000_01_DB.tcl
# the as-flown when unzipped is like
# |2007-077-11:50:02.127500|PPT|Begin|saa-cold-77-10|74566|74566|2|2|33628998|33628998|183.604|14.053|346.321|0|0|10|100.000|Validated|SAA @ 2007-077-11:51:11|00074566002|
# |2007-077-11:53:25.727381|Slew Settled|Begin|saa-cold-77-10|74566|74566|2|2|33628998|33628998|317.522|-30.014|60.266|0|0|10|100.000|Validated|SAA @ 2007-077-11:51:11|00074566002|
# |2007-077-11:55:00.000000|SAA|Begin|||||||||||||||Validated|||
# |2007-077-12:03:10.726820|PPT|End|saa-cold-77-10|74566|74566|2|2|33628998|33628998|317.522|-30.014|60.265|0|0|10|100.000|Validated|SAA @ 2007-077-11:51:11|00074566002|
# The preplanned HTML file is like
# <TR><TD> <FONT size="2">2008-246-03:25:00<TD> <FONT size="2">2008-246-03:37:00<TD> <FONT size="2"><a href = "SDSSJ140642.0+031940_ppt.html">SDSS J140642.0+031940</a><TD> <FONT size="2">211.65111<TD> <FONT size="2">3.3351599<TD> <FONT size="2">287.17465<TD> <FONT size="2">37754<TD> <FONT size="2">1<TD> <FONT size="2">AUTO<TD> <FONT size="2">0x20AB<TD> <FONT size="2">50<TD> <FONT size="2">720
# <TR><TD> <FONT size="2">2008-246-03:37:00<TD> <FONT size="2">2008-246-03:59:00<TD> <FONT size="2"><a href = "saa-cold-246-20_ppt.html">saa-cold-246-20</a><TD> <FONT size="2">313.82042<TD> <FONT size="2">55.000000<TD> <FONT size="2">329.87982<TD> <FONT size="2">74218<TD> <FONT size="2">3<TD> <FONT size="2">AUTO<TD> <FONT size="2">0x0009<TD> <FONT size="2">100<TD> <FONT size="2">1320
# <TR><TD> <FONT size="2">2008-246-03:59:00<TD> <FONT size="2">2008-246-04:29:00<TD> <FONT size="2"><a href = "SGR0501+4516_ppt.html">SGR 0501+4516</a><TD> <FONT size="2">75.278000<TD> <FONT size="2">45.276230<TD> <FONT size="2">81.400000<TD> <FONT size="2">321174<TD> <FONT size="2">11<TD> <FONT size="2">WT<TD> <FONT size="2">0x2227<TD> <FONT size="2">99<TD> <FONT size="2">1800

# The new table format is
# <tr class='header'><td><br />&nbsp;Begin&nbsp;</td><td><br />&nbsp;End&nbsp;</td><td>&nbsp;Target&nbsp;<br />&nbsp;ID&nbsp;</td><td><br />&nbsp;Seg.&nbsp;</td><td><br />&nbsp;Target Name&nbsp;</td><td><br />&nbsp;R.A.&nbsp;</td><td><br />&nbsp;Dec.&nbsp;</td><td><br />&nbsp;Roll&nbsp;</td><td>&nbsp;XRT&nbsp;<br />&nbsp;Mode&nbsp;</td><td>&nbsp;UVOT&nbsp;<br />&nbsp;Mode&nbsp;</td><td><br />&nbsp;FoM&nbsp;</td><td>&nbsp;Time&nbsp;<br />&nbsp;(s)&nbsp;</td></tr>
#		<tr class='norm1'>
#			<td>&nbsp;2009-12-15 00:00:02&nbsp;</td>
#			<td>&nbsp;2009-12-15 00:13:59&nbsp;</td>
#			<td>&nbsp;<a href='https://www.swift.psu.edu/operations/obsSchedule.php?t=30793'>30793</a>&nbsp;</td>
#			<td>&nbsp;<a href='https://www.swift.psu.edu/operations/obsSchedule.php?t=30793&amp;s=103&amp;a=1'>103</a>&nbsp;</td>
#			<td style='text-align:left;'>&nbsp;Mkn501&nbsp;</td>
#			<td>&nbsp;253.473&nbsp;</td>
#			<td>&nbsp;39.792&nbsp;</td>
#			<td>&nbsp;179.670&nbsp;</td>
#			<td>&nbsp;WT&nbsp;</td>
#			<td>&nbsp;0x308f&nbsp;</td>
#			<td>&nbsp;72&nbsp;</td>
#			<td>&nbsp;837&nbsp;</td>


# Old patterns for pointingTable.oldLoad*:
# use datatime.strftime(asflownpattern)
# asflownpattern="/Volumes/Data/Swift/swift-trend/%Y_%m/misc/asflown/AFST_%Y%j*.txt.gz"
# preplannedpattern="/Users/palmer/Documents/networking/downloads/ppst-mail/PPST_%Y%j*.txt*"
# preplannedhtmldir="http://www.swift.psu.edu/operations/PPST/%Y%j/"
# asflownpattern="/Volumes/Data/Swift/swift-trend/%Y_%m/misc/asflown/AFST_%Y%j*.txt.gz"
# preplannedpattern="/Users/palmer/Documents/networking/downloads/ppst-mail/PPST_%Y%j*.txt*"
# preplannedhtmldir="http://www.swift.psu.edu/operations/PPST/%Y%j/"


class pointingTable:
    def __init__(self, ppstfile=None):
        self._daylist = []  # Stored as int(YYYYdoy)
        self._entries = []
        if ppstfile:
            self._explicitPPST = ppstfile
            for f in self._explicitPPST:
                self.loadURLHTML("file://%s" % os.path.realpath(f), 0)
                if verbose:
                    for e in self._entries:
                        print(e)
        else:
            self._explicitPPST = []

    def getPointings(self, trange):
        try:
            tfullrange = (min(trange), max(trange))
        except TypeError:
            tfullrange = (trange, trange)
        if not self._explicitPPST:
            tmin = tfullrange[0] - _maxobslength
            tmax = tfullrange[1] + _maxobslength + datetime.timedelta(days=1,
                                                                      minutes=90)  # Read an extra day for the end of day case
            t = tmin
            while t < tmax:
                try:
                    self.loadDay(t)
                except:
                    if t == tmin:  # First day of range, try stepping back 5 days to see if you hit any useful files
                        if verbose:
                            print("Searching backwards")
                        for lookback in range(-1, -5, -1):
                            try:
                                self.loadDay(t + datetime.timedelta(days=lookback))
                                # print("Found ")
                                break  # If it works, break out
                            except:
                                if verbose:
                                    traceback.print_exc()
                                pass
                t += datetime.timedelta(days=1)
        return self._findEntries(tfullrange)

    def loadDay(self, t):
        day = t.strftime("%Y-%m-%d")
        if verbose:
            print("Loading day %s " % day)
        if day in self._daylist:
            return True  # already loaded
        (daypointlist, asflowncomplete) = self.loadURLHTML(t.strftime(asflownURLpattern), False)
        if not daypointlist or not asflowncomplete:
            if verbose:
                print("%s : %d observations as flown, asflowncomplete = %s" % (day, len(daypointlist), asflowncomplete))
            (ppstdaylist, ppstread) = self.loadURLHTML(t.strftime(preplannedURLpattern), True)
            if not daypointlist:
                daypointlist = ppstdaylist  # No asflown data
            else:
                lastasflown = max([p._tend for p in daypointlist])  # assume there are no gaps in the as-flown
                daypointlist.extend([p for p in ppstdaylist if p._tstart > lastasflown])
        self._entries.extend(daypointlist)
        self._daylist.append(day)

    def loadURLHTML(self, u, preplanned):
        """ Returns (pointlist, complete)
        """
        pointlist = []
        complete = False
        moretocome = False
        if verbose:
            print(u)
        try:
            s = BeautifulSoup.BeautifulSoup(urlopen(u, None, timeout=10, context=unsafe_context), "html.parser")
            table = s.find('table')  # , {"class": "ppst"})  # Even the as-flown table has class ppst
            for row in table.findAll('tr'):
                try:
                    pe = pointingEntry.fromTableRow(row, preplanned)
                    pointlist.append(pe)
                except:
                    text = strip_nbsp(soupCatStrings(row))
                    if verbose:
                        print(text)
                    try:
                        # Crude way of checking whether the AFST is complete.  Better would be to check the time
                        if (text[0:9] == 'There may'
                                or text[0:10] == 'There will'):  # be later observations for this date
                            break  # Do not execute else clause
                    except:
                        pass
            else:
                if len(pointlist) != 0:  # Reached end of table and actually read something
                    # Check if what was read is within 20 min before to 5 minutes after end of day 
                    complete = (pointlist[-1]._tend + datetime.timedelta(minutes=20)).time() < datetime.time(0, 25, 0)
                else:
                    complete = False
        except:
            return (None, False)
        return (pointlist, complete)

    def _findEntries(self, trange):
        try:
            if len(trange) == 1:
                return [p for p in self._entries if p.fContains(trange)]
            else:
                return [p for p in self._entries if p.fOverlaps(trange)]
        except:
            return [p for p in self._entries if p.fContains(trange)]


def soupCatStrings(soupnode, stripit=False):
    """Traverse the tree starting from soupnode, and return the 
    strings of the leaf nodes concatenated together"""
    if isinstance(soupnode, str):
        s = soupnode
    else:
        s = "".join([soupCatStrings(sn, stripit) for sn in soupnode.contents])
    if stripit:
        s.strip()
    return s


class PointingsDatabase:
    """ Database of pointings.  The database is only populated day-by-day as needed"""

    def __init__(self):
        try:
            os.stat(self.Filename)
        except:
            self.CreateFile()

    def ReadDay(self, t, force=False):
        if force or self.CheckDay(t) != 3:
            raise RuntimeError("Not yet implemented")

    def Filename():
        return pdbfile

    Filename = staticmethod(Filename)

    def CheckDay(self, t):
        """ Returns the status of a day in the database:
            0 not loaded
            1 PPST only
            2 PPST and partial AFST
            3 final AFST """
        tstr = t.strftime("%Y-%m-%d")
        conn = sqlite3.connect(pdbfile, detect_types=sqlite3.PARSE_DECLTYPES)
        result = 0
        try:
            daystatus = conn.cursor()
            daystatus.execute('''select asflown_complete,asflow_partial,preplanned from days where date=?''', (tstr,))
            (af, af_partial, preplanned) = daystatus.fetchone()
            if af:
                result = 3
            elif af_partial:
                result = 2
            elif preplanned:
                result = 1
        except:
            pass
        conn.close()
        return result

    def CreateFile(self, delete=False):
        try:
            if delete:
                os.remove(pdbfile)
                raise EOFError("No error")
            else:
                os.stat(pdbfile)
        except:
            conn = sqlite3.connect(pdbfile, detect_types=sqlite3.PARSE_DECLTYPES)
            conn.execute('''create table pointings (%s)''' % pdbfields_text)
            conn.execute('''create table days (%s)''' % pdbdayfields_text)
            conn.commit()
            conn.close()


def hasLength(x):
    try:
        len(x)
        return True
    except:
        return False


def parseNotice(lines):
    """
    Given a source of lines, read the lines as a GRB notice and return (ra_deg,dec_deg,met_swift)
    Where met_swift is adjusted by the Swift UTCF even if the instrument is not Swift.
    If ra and dec are not included in the notice, they are set to None
    """
    ra_deg = None
    dec_deg = None
    strippers = "\xc2\xa0 \t\n"  # some strange characters when Mail.app saves mail
    for l in lines:
        l = swutil.removeNonAscii(l)
        ipntime = ipntimeRE.search(l)
        if ipntime:
            g = ipntime.groupdict()
            ss = float(g['s'])
            h = int(ss / 3600)
            ss -= h * 3600
            m = int(ss / 60)
            ss -= m * 60
            ymdT = "20%s-%s-%sT" % (g['y'], g['m'], g['d'])
            hms = "%02d:%02d:%05.3f" % (h, m, ss)
        else:
            ls = l.split()
            if len(ls) > 1:
                if ls[0] in ["GRB_RA:", "SRC_RA:"]:
                    ra_deg = float(ls[1].split("d")[0].strip(strippers))
                elif ls[0] in ["GRB_DEC:", "SRC_DEC:"]:
                    dec_deg = float(ls[1].split("d")[0].strip(strippers))
                elif ls[0] in ["GRB_DATE:", "DISCOVERY_DATE:"]:
                    ymdT = "20%s-%s-%sT" % tuple(ls[5].strip(strippers).split("/"))
                elif ls[0] in ["GRB_TIME:", "DISCOVERY_TIME:"]:
                    hms = ls[3].strip(strippers + "{}")

                if ls[0] in ["GRB_RA:", "SRC_RA:", "EVENT_RA:"]:
                    ra_deg = float(ls[1].split("d")[0].strip(strippers))
                elif ls[0] in ["GRB_DEC:", "SRC_DEC:", "EVENT_DEC:"]:
                    dec_deg = float(ls[1].split("d")[0].strip(strippers))
                elif ls[0] in ["GRB_DATE:", "EVENT_DATE:", "DISCOVERY_DATE:", "TRIGGER_DATE:"]:
                    ymdT = "20%s-%s-%sT" % tuple(ls[5].strip(strippers).split("/"))
                elif ls[0] in ["GRB_TIME:", "DISCOVERY_TIME:", "EVENT_TIME:", "TRIGGER_TIME:"]:
                    hms = ls[3].strip(strippers + "{}")
    met = swutil.string2met(ymdT + hms, nocomplaint=True, correct=True)
    return (ra_deg, dec_deg, met)


def lineofsight(sat, source):
    """Calculate the
    (elevation_above_Earth_limb_degrees, minimum_elevation_above_sea_level_meters, description)
    of the line of sight from the satellite to the source.
    minimum_elevation_above_sea_level is 0 if LOS intersects Earth, sat.elevation if it LOS is upwards
    """
    zenangle = source.distance(sat.a_ra * r2d, sat.a_dec * r2d).norm
    zenhoriz = ephem.halfpi + math.acos(earthradius_m / (sat.elevation + earthradius_m))
    if zenangle < ephem.halfpi:
        losheight = sat.elevation
        description = "up"
    elif zenangle > zenhoriz:
        losheight = 0.0
        description = "down"
    else:
        losheight = ((sat.elevation + earthradius_m) * math.cos(zenangle - ephem.halfpi)) - earthradius_m
        if losheight > atm_thickness:
            description = "depressed but unattenuated"
        else:
            description = "attenuated by atmosphere at %.0f km" % (losheight / 1000.0,)
    return ((zenhoriz - zenangle) * r2d, losheight, description)


def usage(progname):
    print("Usage: %s MET [MET....]" % progname)
    print("        -v --verbose                    diagnostic information")
    print("           --terse                      produce terse output")
    print("        -o --orbit                      satellite position")
    print("        -s --source sourcename          source visibility.  Sourcename can be '123.4_-56.7' for RA_dec")
    print("        -S --sun --Sun                  Sun visibility")
    print("            --visible                   Only when the source is in the FOV")
    print("        -x --excelvis                   source visibility list in CSV Excelable format (grep vis)")
    # print("        -m --machine                    machine-convenient format with printouts on single lines")
    print("        -p --position \"ra_deg dec_deg\"  position visibility")
    print("        -t --timerange                  use all pointings in the range of times")
    print("           --steptime seconds           cover the range fo times with this interval")
    print("        -c --clipboard                  use times in the current clipboard")
    print("        -f --format '%Y-%m-%dT%H:%M:%S' use given time format.  (Example is default)")
    print("        -P --ppst   ppstfile            use local PPST file instead of getting from web")
    print("        -a --attitude \"ra dec roll\"     manually include an attitude rather than PPST")
    print("           --notice                     process GCN notice piped to stdin to extract time and position")
    # print("           --METonly                    print out nothing but the MET")
    # print("           --UTConly                    print out nothing but the UTC")
    print("        -h --help                       this printout")


# @swutil.timeout(60)
def main(argv=None, debug=None):
    global verbose
    global terse
    global sourcecat
    if argv is None:
        argv = sys.argv
    sources = []
    orbits = False
    ppstfile = []
    timerange = False
    visible_only = False
    timestep = 0
    manualAttitudes = []
    excelvis = False
    sunsource = None
    # debug=open("/tmp/swinfo.debug","a")
    # debug = None

    # if debug :
    #    sys.stdout = swutil.TeeFile(debug, sys.stdout)

    if debug:
        debug.write("%s (%d): %s\n" % (datetime.datetime.now(), os.getpid(), " ".join(argv)))
        vebose = True
        debug.flush()

    format = swutil.fitstimeformat
    try:
        opt, timeargs = getopt.gnu_getopt(argv[1:], "vos:Sxp:tcf:P:a:h",
                                          ["verbose", "terse", "orbit", "source=", "sun", "Sun", "visible",
                                           "excelvis", "position=", "timerange", "steptime=", "clipboard",
                                           "format=", "ppst=", "attitude=", "notice", "help"])
        # print(opt)
        # print(timeargs)
        metvalues = [swutil.string2met(a, correct=True) for a in timeargs]
        for (o, v) in opt:  # option value
            if o in ("-v", "--verbose"):
                verbose = True
                print(" ".join(argv))
            elif o in ("--terse"):
                terse = True
            elif o in ("-o", "--orbit"):
                orbits = orbit()
                # print("Orbit read")
            elif o in ("-s", "--source"):
                loadsourcecat()
                try:
                    # print(v)
                    s = sourcecat.byName(v)
                    # print(s)
                    sources.append(s)
                except:
                    m = radecnameRE.match(v)
                    if m:
                        s = "|||" + v + ("||0||||%(rastring)s|%(decstring)s|" % m.groupdict())
                        sources.append(source(s))
                    else:
                        print("Source %s not found in catalog" % v)
                        return
            elif o in ("-S", "--sun", "--Sun"):
                s = source("Sun")
                sources.append(s)
            elif "--visible".startswith(o):
                visible_only = True
            elif o in ("-t", "--timerange"):
                timerange = True
            elif o in ("--steptime",):
                timestep = datetime.timedelta(seconds=float(v))
            elif o in ("-c", "--clipboard"):
                for ttry in os.popen("/usr/bin/pbpaste | tr '\r' '\n'"):
                    try:
                        t = swutil.string2met(ttry, nocomplaint=True, correct=True)
                        if t:
                            metvalues.append(t)
                    except:
                        pass
            elif o in ("-f", "--format"):
                format = v
            elif o in ("-p", "--position"):
                try:
                    ra, dec = [ephem.degrees(s) for s in v.translate(split_translate).split()]
                    string = "|||%s||%d||||%f|%f|" % (v, len(sources), ra / ephem.degree, dec / ephem.degree)
                    newsource = source(string)
                    sources.append(newsource)
                except:
                    print("Could not add source at position " + v)
                    return
            elif o in ("-P", "--ppst"):
                ppstfile.append(v)
            elif o in ("-a", "-attitude"):
                vsplit = v.translate(split_translate).split()
                try:
                    ra = float(vsplit[0])
                    dec = float(vsplit[1])
                    roll = float(vsplit[2])
                except:
                    print("Could not get ra,dec,roll out of argument '%s'" % (v,))
                    usage(argv[0])
                    return
                ztime = datetime.datetime(2000, 1, 1, 0, 0, 0)
                manualAttitudes.append(
                    pointingEntry(ztime, ztime, ztime, ra, dec, roll, 0, 0, "(from command line)", False))
            elif o in ("--notice",):
                (ra, dec, met) = parseNotice(sys.stdin)
                metvalues.append(met)
                if ra != None and dec != None:
                    string = "|||%.2f_%.2f||%d||||%f|%f|" % (ra, dec, len(sources), ra, dec)
                    newsource = source(string)
                    sources.append(newsource)
            elif o in ("-x", "--excelvis"):
                excelvis = True
            elif o in ("-h", "--help"):
                usage(argv[0])
                return

        if debug: debug.write("Arguments processed\n")

        thePointingTable = pointingTable(ppstfile)

        if debug: debug.write("Pointing table received \n")

        if timerange:
            pytimes = (swutil.met2datetime(min(metvalues), correct=True),
                       swutil.met2datetime(max(metvalues), correct=True))
            print("Time range = %s - %s " % pytimes)
            # print(thePointingTable.getPointings(pytimes))
            # Get the METs of the middles of the observations
            metvalues = [swutil.datetime2met(p._tstart +
                                             datetime.timedelta(seconds=
                                                                (p._tend - p._tstart).total_seconds() / 2.0),
                                             correct=True)
                         for p in thePointingTable.getPointings(pytimes)]

        if timestep:
            pytimes = (swutil.met2datetime(min(metvalues), correct=True),
                       swutil.met2datetime(max(metvalues), correct=True))
            print("Time range = %s - %s by %f seconds" % (pytimes + (timestep.total_seconds(),)))

            metvalues = [swutil.datetime2met(pytimes[0] + i * timestep, correct=True)
                         for i in range(1 + int(math.ceil((pytimes[1] - pytimes[0]).total_seconds()
                                                          / timestep.total_seconds())))]
            # Get the METs of the middles of the observations

        if excelvis:
            print("vis,Source,slewstart,slewend,obsend,exposure,elevstart,elevend")
            # Excel visibility output requires orbit information
            # to calculate elevation at start of observation=end of slew, end of observation
            if not orbits:
                orbits = orbit()

        if verbose and sources:
            print(sources)

        if debug:
            debug.write("About to loop through metvalues: %s\n" % (metvalues,))
            debug.flush()

        for t in metvalues:
            pointprint = StringIO()
            visible = len(sources) == 0  # even if visible_only, print if no sources
            u = utcf(t, verbose)
            udelta = datetime.timedelta(seconds=u)
            if hasLength(u):
                ulist = u
                tlist = t
                pytime = [swutil.met2datetime(t_, correct=True) for t_ in t]
                pystarttime = pytime[0]
            else:
                ulist = [u]
                tlist = [t]
                pytime = swutil.met2datetime(t, correct=True)
                pystarttime = pytime
            for i in range(len(ulist)):
                if not terse:
                    print("Time(MET + UTCF -> UT):   ", file=pointprint, end='')
                print("%.3f + %.3f -> %s" % (tlist[i] + 5e-4, ulist[i],
                                             swutil.met2fitsString(tlist[i], milliseconds=True, correct=True,
                                                                   format=format)), file=pointprint)
                if not terse:
                    print(swutil.met2fitsString(tlist[i], milliseconds=True, correct=True,
                                                format="YYMMDD_SOD:                   %y%m%d_%q       DOY=%j"),
                          file=pointprint)
            if manualAttitudes:
                pointings = manualAttitudes
            else:
                pointings = thePointingTable.getPointings(pytime)
            if orbits:
                if debug:
                    debug.write("About to get the satellite position for %s\n" % (pystarttime,))
                    debug.flush()
                sat = orbits.getSatellite(pystarttime)
                if debug:
                    debug.write("Got the satellite position\n")
                    debug.flush()
                if terse:
                    print("Satellite zenith RA=%.2f, dec=%.2f, lat=%.2f N, lon=%.2f E" % (
                        sat.a_ra * r2d, sat.a_dec * r2d, sat.sublat * r2d, sat.sublong * r2d))
                else:
                    print("Swift Zenith(RA,dec):         %.2f, %.2f" % (sat.a_ra * r2d, sat.a_dec * r2d))
                    print("Swift Location(lon,lat,alt):  %.2f E, %.2f N, %.0f km" % (
                        sat.sublong * r2d, sat.sublat * r2d, sat.elevation / 1000))
            for p in pointings:
                if not excelvis:
                    if terse:
                        print(p, file=pointprint)
                    else:
                        print("Obs Sequence Number:          %08d%03d" % (p._obs, p._segment), file=pointprint)
                        print("Obs Target Name:              %s%s" % (
                            p._sourcename, (" (planned)" if p._planned else "")), file=pointprint)
                        print("Obs Date and Times:           %s - %s" % (
                            (p._tstart + udelta).strftime("%Y-%m-%d   %H:%M:%S"),
                            (p._tend + udelta).strftime("%H:%M:%S")),
                              file=pointprint)
                        print("Obs Pointing(ra,dec,roll):    %.3f, %.3f, %.2f" % (p._ra, p._dec, p._roll),
                              file=pointprint)
                for s in sources:
                    try:
                        s.compute(pystarttime)
                        theta, phi = s.thetangle_phi(p._ra, p._dec, p._roll)
                        # print("exp,cosangle = ",s.exposure(p._ra, p._dec, p._roll) )
                        exp, cosangle = s.exposure(p._ra, p._dec, p._roll)
                        visible = visible or (exp > 0)
                        if orbits:
                            (elev, losheight, description) = lineofsight(sat, s)
                            relhoriz = " (%s)" % (description,)
                        else:
                            relhoriz = ""
                        if terse:
                            print("%s: (theta,phi) = (%.4f, %.4f); exposure = %.0f cm^2%s" % (
                                s.name, theta, ephem.degrees(phi).znorm, exp * cosangle, relhoriz), file=pointprint)
                        elif excelvis:
                            (elevstart, losheightstart, descriptionstart) = lineofsight(
                                orbits.getSatellite(p._tslewend), s)
                            (elevend, losheightend, descriptionend) = lineofsight(orbits.getSatellite(p._tend), s)
                            print("vis,%s,%s,%s,%s,%.0f,%.1f,%.1f" % (s.name,
                                                                      p._tstart.strftime("%Y-%m-%d %H:%M:%S"),
                                                                      p._tslewend.strftime("%Y-%m-%d %H:%M:%S"),
                                                                      p._tend.strftime("%Y-%m-%d %H:%M:%S"),
                                                                      exp * cosangle,
                                                                      elevstart,
                                                                      elevend), file=pointprint)
                        else:
                            print("%s_imageloc (boresight_dist, angle): (%.0f, %.0f)" % (
                                s.name, theta * r2d, ephem.degrees(phi).znorm * r2d), file=pointprint)
                            print("%s_exposure (cm^2 cos adjusted): %.0f" % (s.name, exp * cosangle), file=pointprint)
                            if len(relhoriz) > 3:
                                print("%s_altitude: %.0f (%s)" % (s.name, elev, description), file=pointprint)
                    except:
                        traceback.print_exc()
                        # traceback.print_tb(sys.exc_info()[2])
                        print(s.catnum + "ran into trouble in calculations")
                        return
            if visible or not visible_only:
                sys.stdout.write(pointprint.getvalue())
    except:
        usage(argv[0])
        traceback.print_exc()
        if debug: traceback.print_exc(None, debug)
        # traceback.print_tb(sys.exc_info()[2])
    if debug:
        debug.write("%s: swinfo main() completed\n" % (datetime.datetime.now(),))
        debug.flush()


def checkParse():
    """ Go through all the parsed mail messages to see which ones weren't parsed properly"""
    for f in glob.glob('/tmp/latesttrigger.*.processed'):
        parsed = parseNotice(open(f).readlines())
        if len(parsed) != 3 or parsed[2] < 3e8:
            print("%s %s" % (f, parsed))


if __name__ == "__main__":
    # for i in range(len(sys.argv)) :
    #    print("%d: %s" % (i, sys.argv[i]))
    # Turn on debugging
    debug = None
    debug = open("/tmp/swinfo.debug", "a")
    # if debug:
    #     swutil.dumponsignal(fname="/tmp/swinfo.debug")
    #     main = swutil.TeeStdoutDecorator(main, debug)
    if len(sys.argv) > 1:
        main(sys.argv, debug=debug)
    else:
        main([l for l in os.popen("/usr/bin/pbpaste | tr '\r' '\n'")])
