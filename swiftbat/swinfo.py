#! /usr/bin/env python

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
import functools
from pathlib import Path
import re
import datetime
import shutil
import glob
import gzip
import traceback
import getopt
import time
from . import swutil
from .clockinfo import utcf
from . import sftime, sfts, loadsfephem
from .sfmisc import _asutc
import astropy.units as u
from astropy.io import fits
from astropy import coordinates
from astropy.coordinates import ICRS, SkyCoord, Angle as apAngle
import numpy as np
from io import StringIO
from functools import lru_cache as __lru_cache
import swifttools.swift_too as swto
from urllib.request import urlopen, quote
from swiftbat import generaldir
import skyfield.api as sf_api
from skyfield.trigonometry import position_angle_of as sf_position_angle_of
from skyfield.positionlib import ICRF as sf_ICRF
import sqlite3  # Good sqlite3 tutorial at http://www.w3schools.com/sql/default.asp
from typing import List, Tuple


split_translate = str.maketrans("][, \t;", "      ")
# Skyfield interfaces
_sf_load = sf_api.Loader("~/.skyfield", verbose=False)
_sf_ts = _sf_load.timescale(builtin=True)


@__lru_cache(0)
def loadsfephem():
    return _sf_load("de421.bsp")


_sf_timescale = sf_api.load.timescale(builtin=True)


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

execdir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))

basecatalog = os.path.join(execdir, "catalog")
fitscatalog = os.path.join(execdir, "recent_bcttb.fits.gz")
# FIXME should be handled by dotswift
catalogdir = "/opt/data/Swift/analysis/sourcelist"
newcatalogfile = os.path.join(catalogdir, "newcatalog")
cataloglist = [
    os.path.join(catalogdir, "trumpcatalog"),
    os.path.join(catalogdir, "grbcatalog"),
    os.path.join(catalogdir, "catalog"),
    newcatalogfile,
]


# define machineReadable=False

ydhms = "%Y-%j-%H:%M:%S"

# TLEpattern = ["ftp://heasarc.gsfc.nasa.gov/swift/data/obs/%Y_%m/",".*","auxil","SWIFT_TLE_ARCHIVE.*.gz"]
TLEpattern = [
    "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/%Y_%m/",
    ".*",
    "auxil",
    "SWIFT_TLE_ARCHIVE.*.gz",
]

tlefile = "/tmp/latest_swift_tle.gz"
tlebackup = os.path.expanduser("~/.swift/recent_swift_tle.gz")

radecnameRE = re.compile(r"""^(?P<rastring>[0-9.]+)_+(?P<decstring>([-+]*[0-9.]+))$""")

ipntimeRE = re.compile(
    r"""'(?P<mission>[^']+)'\s*'(?P<d>[ 0-9]+)/(?P<m>[ 0-9]+)/(?P<y>[ 0-9]+)'\s+(?P<s>[0-9.]+)"""
)


verbose = False  # Verbose controls diagnositc output
terse = False  # Terse controls the format of ordinary output


# Old fashioned
def simbadnames(query):
    """Given a source name, or other SIMBAD query, generates a list of identifier matches
    See http://simbad.u-strasbg.fr/simbad/sim-fscript
    """
    url_ = urlopen(
        """http://simbad.u-strasbg.fr/simbad/sim-script?submit=submit+script&script=format+object+%%22+%%25IDLIST%%22%%0D%%0Aquery++%s"""
        % (quote(query),),
        None,
        60,
    )
    names = []
    while url_:
        l = url_.readline().decode()
        if re.match("""^::error::*""", l):
            raise ValueError(query)
        if re.match("""^::data::*""", l):
            break
    for l in url_.readlines():
        s = l.decode().strip().split()
        if len(s) > 0:
            if s[0] == "NAME":
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

        if "RA" in table.keys():
            ra_string = "RA"
            dec_string = "DEC"
        else:
            ra_string = "ra"
            dec_string = "dec"

        co = coordinates.SkyCoord(
            table[ra_string][0], table[dec_string][0], unit=(u.hour, u.deg), frame="fk5"
        )
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
        allTLE = gzip.open(
            tlefile, "rt" if sys.version_info.major == 3 else "r"
        ).readlines()
        # Some TLE files seem to have additional blank lines
        # Actually, what they have is \r\n in some cases and \n in others
        # Nevertheless, put this in in case they change the format
        allTLE = [l.strip() for l in allTLE if len(l.strip()) == 69]
        nTLE = len(allTLE) // 2
        self._tleByTime = [
            (
                datetime.datetime(
                    2000 + int(allTLE[2 * i][18:20]), 1, 1, 0, 0, 0, tzinfo=datetime.UTC
                )
                + datetime.timedelta(float(allTLE[2 * i][20:32]) - 1),
                ("Swift", allTLE[2 * i], allTLE[2 * i + 1]),
            )
            for i in range(nTLE)
        ]
        self.pickTLE(datetime.datetime.now(datetime.UTC))

    def pickTLE(self, t):
        t = _asutc(t)
        if self._tleByTime[-1][0] < t:
            self._tle = self._tleByTime[-1][1]
            self._tletimes = [
                self._tleByTime[-1][0],
                datetime.datetime(datetime.MAXYEAR, 1, 1, tzinfo=datetime.timezone.utc),
            ]
        else:
            for w in range(len(self._tleByTime) - 1, -1, -1):
                if self._tleByTime[w][0] < t or w == 0:
                    self._tle = self._tleByTime[w][1]
                    self._tletimes = [self._tleByTime[w][0], self._tleByTime[w + 1][0]]
                    break

    def getSatellite(self, t) -> Tuple[sf_api.EarthSatellite, sf_api.Time]:
        """Produces skyfield.EarthSatellite and its location for Swift at the given time

        Args:
            t (_type_): _description_

        Returns:
            _type_: _description_
        """
        global verbose
        t = _asutc(t)
        # print(t, self._tletimes)
        if t < self._tletimes[0] or self._tletimes[1] < t:
            if verbose:
                print("Picking TLE for time %s" % (t,))
            self.pickTLE(t)
        if verbose:
            print(
                "Picked TLE for %s < %s < %s"
                % (self._tletimes[0], t, self._tletimes[1])
            )
        sat = sf_api.EarthSatellite(self._tle[1], self._tle[2], name=self._tle[0])
        sft = sftime(t)
        return sat, sat.at(sft)

    def satname(self):
        return self._tleByTime[0][1][0]

    def updateTLE(self, maxage_days=10):
        global verbose
        try:
            if not Path(tlefile).exists():
                if Path(tlebackup).exists():
                    shutil.copyfile(tlebackup, tlefile)
                    shutil.copystat(tlebackup, tlefile)
            checksecs = time.mktime(
                (
                    datetime.datetime.now(datetime.UTC)
                    + datetime.timedelta(days=-maxage_days)
                ).timetuple()
            )
            if os.stat(tlefile).st_mtime > checksecs:
                return  # The TLE file exists and is less than maxage
        except:
            pass
        tlematch = re.compile(TLEpattern[-1])
        for url_month in (
            datetime.datetime.now(datetime.UTC).strftime(TLEpattern[0]),
            (datetime.datetime.now(datetime.UTC) - datetime.timedelta(30)).strftime(
                TLEpattern[0]
            ),
        ):
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
                                # Cache the most recent TLE file
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
    theta = apAngle(theta, u.rad)
    phi = apAngle(phi, u.rad)

    if np.cos(theta) < 0:
        return (0.0, np.cos(theta).value)

    if theta > apAngle(90, u.deg):
        return (0.0, 0.0)

    # BAT dimensions
    detL = 286 * 4.2e-3  # Amynote: each det element is has length 4.2e-3 m
    detW = 173 * 4.2e-3  # Amynote: each det element is has length 4.2e-3 m
    maskL = 487 * 5.0e-3  # Using the 5 mm mask cells, true size is 95.9" x 47.8"
    maskW = 243 * 5.0e-3
    efl = 1.00
    # Calculate dx, but reverse it if necessary to put projected detector to left of mask center
    dx = (maskL - detL) / 2 - efl * np.tan(theta) * abs(np.cos(phi))
    dy = (maskW - detW) / 2 + efl * np.tan(theta) * np.sin(phi)
    # Boundaries of the detector as clipped to rectangle of mask
    x1 = max(0, -dx)
    x2 = min(detL, maskL - dx)
    deltaX = x2 - x1
    y1 = max(0, -dy)
    y2 = min(detW, maskW - dy)
    deltaY = y2 - y1
    # Now adjust for the cut corners, which first attacks the upper left detector corner
    # as described by the (if-necessary-reversed) coord. system
    xint = (
        (y1 + dy) - maskW / 2 - (dx + x1)
    )  # line of corner clip extends through (x1 + xint, y1)
    if deltaX < 0 or deltaY < 0:
        area = 0
    elif xint <= -deltaY:  # no clipping
        area = deltaX * deltaY
    elif xint <= 0:  # upper left corner clipped along left edge
        if deltaY <= deltaX - xint:  # clip from left edge, to top edge
            area = deltaX * deltaY - ((deltaY + xint) ** 2) / 2
        else:  # clip from left edge to right edge
            area = deltaX * -xint + (deltaX**2) / 2
    elif xint <= deltaX:  # clipped at bottom edge
        if xint <= deltaX - deltaY:  # clipped bottom and top
            area = (deltaX - xint) * deltaY - (deltaY**2) / 2
        else:  # clipped bottom and right
            area = ((deltaX - xint) ** 2) / 2
    else:
        area = 0
    # if you want to see what the corners do: area = max(0,deltaX * deltaY) - area
    # multiply by 1e4 for cm^2, 1/2 for open area
    area = area * 1e4 * u.cm**2 / 2

    return (area.value, np.cos(theta).value)


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
    # bit7 is left or right half
    bit6, bit5_0 = np.divmod(det_in_mod, np.int16(64))
    bit6_3, bit2_0 = np.divmod(
        det_in_mod, np.int16(8)
    )  # bits 6-3 determine the imod, bits 2-0 jmod
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
    # 16 detectors + 2-space gap for second column
    iblock = imod + np.int16(18) * bit10
    jmodrow = (np.array([6, 7, 4, 5, 2, 3, 0, 1]) * 11).astype(
        np.int16
    )  # 8 detectors and 3-space gap for each row of modules
    jblock = jmod + jmodrow[bit9_7]

    # blocks are arranged:
    #       0  1  2  3  4  5  6  7
    #       8  9  10 11 12 13 14 15
    # and blocks 8-15 are rotated 180 degrees
    bit15, bit14_11 = np.divmod(block, np.int16(8))
    #                      8-15        0-7
    y = np.where(bit15, np.int16(85) - jblock - 1, np.int16(88) + jblock)
    x = np.where(bit15, np.int16(34) - iblock - 1, iblock) + np.int16(36) * bit14_11
    if scalar:
        x = np.int16(x)
        y = np.int16(y)
    return x, y


@functools.lru_cache(maxsize=0)
def xy2detidmap():
    """
    Produce a detector map filled with detector IDs

    -1 for unpopulated detector locations

    :return: detids[y, x]
    :rtype: uint16, shape = (173, 286)
    """
    detids = np.arange(0, 2**15)
    x, y = detid2xy(detids)
    result = np.full((y.max() + 1, x.max() + 1), np.int16(-1))
    result[y, x] = detids
    return result


def xy2detid(x, y):
    dmap = xy2detidmap()
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
            sourcecat = sourcelist(
                cataloglist[0:3] + [fitscatalog], verbose=verbose
            )  # Exclude newcatalog
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

        :param initstring: String from source table ('|' delimited), or object recognized by Simbad
        :param kwargs: {ra:ra_deg, dec:dec_deg, <name:'a name'>, <catnum:catnum>}
        """
        # This is where the Sun was implemented
        # if initstring in ephem.__dict__:
        #     self.needs_computing = True
        #     self.computable = ephem.__dict__[initstring]()
        # elif...
        if "ra" in kwargs and "dec" in kwargs:
            if initstring:
                raise RuntimeError("Give ra=,dec= or an initstring but not both")
            self.set_radec(kwargs["ra"], kwargs["dec"])
            self.name = kwargs.get("name", "unnamed")
            self.catnum = kwargs.get("catnum", 0)
        else:
            self.needs_computing = False
            if "|" in initstring:
                s = initstring.split("|")
                self.name = s[3].strip().replace(" ", "_")
                self.catnum = int(s[5])
                self.set_radec(float(s[9]), float(s[10]))
            else:
                try:
                    ra, dec = simbadlocation(initstring)
                    self.set_radec(ra, dec)
                    self.name = initstring
                    self.catnum = -1
                except:
                    raise NameError(
                        f"Source name {initstring} not recognized by Simbad"
                    )

    def set_radec(self, ra_deg, dec_deg):
        self.skyloc = SkyCoord(
            ICRS(ra=apAngle(ra_deg, u.deg), dec=apAngle(dec_deg, u.deg))
        )

    @property
    def ra_deg(self):
        return self.skyloc.ra.degree

    @property
    def dec_deg(self):
        return self.skyloc.dec.degree

    @classmethod
    def source(cls, ra, dec, name="anonymous", catnum=-1) -> str:
        """The source as a '|'-delimted string as used in ASCII catalog tables

        Args:
            ra (float): degrees
            dec (float): degrees
            name (str, optional): _description_. Defaults to "anonymous".
            catnum (int, optional): _description_. Defaults to -1.

        Returns:
            str: _description_
        """
        return cls(
            f"|||{name}||{catnum}||||{apAngle(ra, u.deg).degree}|{apAngle(dec, u.deg).degree}"
        )

    def sf_position(self) -> sf_ICRF:
        """Position of source as a skyfield.Position"""
        return sf_api.position_of_radec(
            ra_hours=self.ra_deg / 15, dec_degrees=self.dec_deg
        )

    def exposure(self, ra, dec, roll):
        # returns (projected_area, cos(theta))
        (thetangle, phi) = self.thetangle_phi(ra, dec, roll)

        return batExposure(thetangle, phi)

    def thetangle_phi(
        self, ra: float, dec: float, roll: float
    ) -> Tuple[(apAngle, apAngle)]:
        """_Source position in instrument FOV given instrument pointing direction
            returns (thetangle,phi) where thetangle is the angular distance from
            the boresight and phi is the angle from the phi=0 axis of BAT.

            theta = tan(thetangle) gives the theta we use, which is the projected distance to a flat plane,
            but it is not useful for thetangle > 90 degrees

        Args:
            ra (float): ra of spacecraft boresight in degrees
            dec (float): dec of spacecraft boresight in degrees
            roll (float): roll is 'roll left' (CCW as seen from behind looking out along boresight) in degrees
        Returns:
            (theta:degree, phi:apAngle): location of source in spacecraft spherical coordinates
        """
        # Use astropy coordinates
        boreloc = SkyCoord(ICRS(ra=apAngle(ra, u.deg), dec=apAngle(dec, u.deg)))
        thetangle = boreloc.separation(self.skyloc).to(u.deg)
        # posang is CCW, phi is CCW from +Y so (posang - roll - 90deg) gives phi
        phi = (
            (
                coordinates.position_angle(
                    boreloc.ra, boreloc.dec, self.skyloc.ra, self.skyloc.dec
                )
                - apAngle(roll, u.deg)
                - apAngle(90, u.deg)
            )
            .wrap_at(180 * u.deg)
            .to(u.deg)
        )
        return (thetangle, phi)

    def compute(self, t):
        """Default implementation does nothing

        Args:
            t (time): _description_
        """
        pass

    def __str__(self):
        return "|||%s||%d||||%f|%f|" % (
            self.name,
            self.catnum,
            self.skyloc.ra.degree,
            self.skyloc.dec.degree,
        )


class sunsource(source):
    def __init__(self, t=None, bodyname="Sun"):
        self.body = loadsfephem()[bodyname]
        self.earth = loadsfephem()["Earth"]
        skyloc = self._skyloc_at(t)
        super().__init__(ra=skyloc.ra.degree, dec=skyloc.dec.degree, name=bodyname)

    def _skyloc_at(self, t=None) -> ICRS:
        # Cached so it only downloads the DE422 once and only reads it once per run
        if t is None:
            t = datetime.datetime.now(tz=datetime.timezone.utc)
        t = sftime(t)
        app = self.earth.at(t).observe(self.body).apparent()
        ra, dec, _ = app.radec()
        skyloc = SkyCoord(
            ICRS(ra=apAngle(ra._degrees, u.deg), dec=apAngle(dec._degrees, u.deg))
        )
        return skyloc

    def compute(self, t):
        self.skyloc = self._skyloc_at(t)


# This is a truncated version of catalog.  It should become the parent of catalog someday (FIXME 2008-09-17)
class sourcelist:
    def __init__(self, filelist=None, verbose=False):
        if filelist is None:
            filelist = [
                os.path.join(catalogdir, "catalog"),
                os.path.join(catalogdir, "grbcatalog"),
            ]
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
                aSource = source.source(
                    ra=row.field("RA_OBJ"),
                    dec=row.field("DEC_OBJ"),
                    name=row.field("NAME"),
                )
                if verbose:
                    print(aSource)
                # print(line,aSource.name)
                self.allsources[self.canonName(aSource.name)] = aSource
        else:
            for line in open(file).readlines():
                line = line.strip()
                if (
                    "+------+------------+" in line
                    or "|   ROW|        TIME|           NAME|" in line
                    or line.startswith("#")
                    or line == ""
                ):
                    continue
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
            ss = float(g["s"])
            h = int(ss / 3600)
            ss -= h * 3600
            m = int(ss / 60)
            ss -= m * 60
            ymdT = "20%s-%s-%sT" % (g["y"], g["m"], g["d"])
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
                elif ls[0] in [
                    "GRB_DATE:",
                    "EVENT_DATE:",
                    "DISCOVERY_DATE:",
                    "TRIGGER_DATE:",
                ]:
                    ymdT = "20%s-%s-%sT" % tuple(ls[5].strip(strippers).split("/"))
                elif ls[0] in [
                    "GRB_TIME:",
                    "DISCOVERY_TIME:",
                    "EVENT_TIME:",
                    "TRIGGER_TIME:",
                ]:
                    hms = ls[3].strip(strippers + "{}")
    met = swutil.string2met(ymdT + hms, nocomplaint=True, correct=True)
    return (ra_deg, dec_deg, met)


def lineofsight(satpos, source):
    """Calculate the
    (elevation_above_Earth_limb_degrees, losheight_m, description)
    of the line of sight from the satellite position to the source.
    minimum_elevation_above_sea_level is 0 if LOS intersects Earth, sat.elevation if LOS is upwards
    """
    atm_thickness_m = 100e3  # How deep in the atmosphere you should report attenuation
    # Using Skyfield positions for separation, thereafter raw radians
    zenangle_rad = satpos.separation_from(source.sf_position()).radians
    _, _, satradius = satpos.radec()
    satradius_m = satradius.m
    earthrad_m = u.earthRad.to(u.m)
    # Hard earth horizon (equatorial radius)
    zenhoriz_rad = np.pi / 2 + np.arccos(earthrad_m / satradius_m)
    if zenangle_rad < np.pi / 2:
        losheight_m = satradius_m - earthrad_m
        description = "up"
    elif zenangle_rad > zenhoriz_rad:
        losheight_m = 0.0
        description = "down"
    else:
        losheight_m = (satradius_m * np.cos(zenangle_rad - np.pi / 2)) - earthrad_m
        if losheight_m > atm_thickness_m:
            description = "depressed but unattenuated"
        else:
            description = "attenuated by atmosphere at %.0f km" % (
                losheight_m / 1000.0,
            )
    return (np.rad2deg(zenhoriz_rad - zenangle_rad), losheight_m, description)


class PointingEntry:
    def __init__(
        self, tstart, tslewend, tend, ra, dec, roll, obs, segment, sourcename, planned
    ):
        self._tstart = tstart
        self._tslewend = tslewend
        if tend == None:  # How to handle the stubs
            self._tend = tstart + datetime.timedelta(seconds=90 * 60)
        else:
            self._tend = tend
        self._ra = ra
        self._dec = dec
        self._roll = roll
        self._obs = obs
        self._segment = segment
        self._sourcename = sourcename
        self._planned = planned

    @classmethod
    def listfromquery(cls, query: swto.ObsQuery) -> List["PointingEntry"]:
        result = []
        for entry in query:
            # query.to_utctime()
            pointent = cls(
                entry.begin,
                entry.begin + entry.slewtime,
                entry.end,
                entry.ra,
                entry.dec,
                entry.roll,
                int(entry.obsid[:-3]),
                entry.seg,
                entry.targname,
                False,
            )
            result.append(pointent)
        return result

    @classmethod
    def listfortimes(cls, times: list[datetime.datetime], allinrange=False):
        if allinrange:
            begin = min(times)
            end = max(times)
            queries = [swto.ObsQuery(begin=begin, end=end)]
        else:
            queries = [swto.ObsQuery(t) for t in times]

        result = []
        for query in queries:
            if query.submit():
                result.extend(cls.listfromquery(query))
        return result

    def __str__(self):
        if self._planned:
            if self._planned == True:
                plannedstring = "(Preplanned)"
            else:
                plannedstring = "(planned %s)" % self._planned
        else:
            plannedstring = ""
        return "%08d%03d  %s - %s - %s [%7.3f, %7.3f, %7.3f] %s %s" % (
            self._obs,
            self._segment,
            self._tstart.strftime("%Y-%j-%H:%M:%S"),
            self._tslewend.strftime("%H:%M:%S"),
            self._tend.strftime("%H:%M:%S"),
            self._ra,
            self._dec,
            self._roll,
            self._sourcename,
            plannedstring,
        )


def usage(progname):
    print("Usage: %s MET [MET....]" % progname)
    print("        -v --verbose                    diagnostic information")
    print("           --terse                      produce terse output")
    print("        -o --orbit                      satellite position")
    print(
        "        -s --source sourcename          source visibility.  Sourcename can be '123.4_-56.7' for RA_dec"
    )
    print("        -S --sun --Sun                  Sun visibility")
    print("            --visible                   Only when the source is in the FOV")
    print(
        "        -x --excelvis                   source visibility list in CSV Excelable format (grep vis)"
    )
    # print("        -m --machine                    machine-convenient format with printouts on single lines")
    print('        -p --position "ra_deg dec_deg"  position visibility')
    print(
        "        -t --timerange                  use all pointings in the range of times"
    )
    print(
        "           --steptime seconds           cover the range fo times with this interval"
    )
    print("        -c --clipboard                  use times in the current clipboard")
    print(
        "        -f --format '%Y-%m-%dT%H:%M:%S' use given time format.  (Example is default)"
    )
    print(
        "        -P --ppst   ppstfile            use local PPST file instead of getting from web"
    )
    print(
        '        -a --attitude "ra dec roll"     manually include an attitude rather than PPST'
    )
    print(
        "           --notice                     process GCN notice piped to stdin to extract time and position"
    )
    # print("           --METonly                    print out nothing but the MET")
    # print("           --UTConly                    print out nothing but the UTC")
    print("        -h --help                       this printout")


# @swutil.timeout(60)
def swinfo_main(argv=None, debug=None):
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
    # debug=open("/tmp/swinfo.debug","a")
    # debug = None

    # if debug :
    #    sys.stdout = swutil.TeeFile(debug, sys.stdout)

    if debug:
        debug.write(
            "%s (%d): %s\n" % (datetime.datetime.now(), os.getpid(), " ".join(argv))
        )
        vebose = True
        debug.flush()

    format = swutil.fitstimeformat
    try:
        opt, timeargs = getopt.gnu_getopt(
            argv[1:],
            "vos:Sxp:tcf:P:a:h",
            [
                "verbose",
                "terse",
                "orbit",
                "source=",
                "sun",
                "Sun",
                "visible",
                "excelvis",
                "position=",
                "timerange",
                "steptime=",
                "clipboard",
                "format=",
                "ppst=",
                "attitude=",
                "notice",
                "help",
            ],
        )
        # print(opt)
        # print(timeargs)
        metvalues = [swutil.string2met(a, correct=True) for a in timeargs]
        for o, v in opt:  # option value
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
                        s = (
                            "|||"
                            + v
                            + ("||0||||%(rastring)s|%(decstring)s|" % m.groupdict())
                        )
                        sources.append(source(s))
                    else:
                        print("Source %s not found in catalog" % v)
                        return
            elif o in ("-S", "--sun", "--Sun"):
                s = sunsource()
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
                    ra, dec = [float(s) for s in v.translate(split_translate).split()]
                    string = "|||%s||%d||||%f|%f|" % (v, len(sources), ra, dec)
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
                    PointingEntry(
                        ztime,
                        ztime,
                        ztime,
                        ra,
                        dec,
                        roll,
                        0,
                        0,
                        "(from command line)",
                        False,
                    )
                )
            elif o in ("--notice",):
                (ra, dec, met) = parseNotice(sys.stdin)
                metvalues.append(met)
                if ra != None and dec != None:
                    string = "|||%.2f_%.2f||%d||||%f|%f|" % (
                        ra,
                        dec,
                        len(sources),
                        ra,
                        dec,
                    )
                    newsource = source(string)
                    sources.append(newsource)
            elif o in ("-x", "--excelvis"):
                excelvis = True
            elif o in ("-h", "--help"):
                usage(argv[0])
                return

        if debug:
            debug.write("Arguments processed\n")

        # thePointingTable = pointingTable(ppstfile)

        if debug:
            debug.write("Pointing table received \n")

        if timerange:
            pytimes = (
                swutil.met2datetime(min(metvalues), correct=True),
                swutil.met2datetime(max(metvalues), correct=True),
            )
            print("Time range = %s - %s " % pytimes)
            # print(thePointingTable.getPointings(pytimes))
            # Get the METs of the middles of the observations
            pointings = PointingEntry.listfortimes(pytimes, allinrange=True)
            metvalues = [
                swutil.datetime2met(
                    p._tstart
                    + datetime.timedelta(
                        seconds=(p._tend - p._tstart).total_seconds() / 2.0
                    ),
                    correct=True,
                )
                for p in pointings
            ]

        if timestep:
            pytimes = (
                swutil.met2datetime(min(metvalues), correct=True),
                swutil.met2datetime(max(metvalues), correct=True),
            )
            print(
                "Time range = %s - %s by %f seconds"
                % (pytimes + (timestep.total_seconds(),))
            )

            metvalues = [
                swutil.datetime2met(pytimes[0] + i * timestep, correct=True)
                for i in range(
                    1
                    + int(
                        np.ceil(
                            (pytimes[1] - pytimes[0]).total_seconds()
                            / timestep.total_seconds()
                        )
                    )
                )
            ]
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
            # even if visible_only, print if no sources
            visible = len(sources) == 0
            utcf_ = utcf(t, verbose)
            udelta = datetime.timedelta(seconds=utcf_)
            if hasLength(utcf_):
                ulist = utcf_
                tlist = t
                pytime = [swutil.met2datetime(t_, correct=True) for t_ in t]
                pystarttime = pytime[0]
            else:
                ulist = [utcf_]
                tlist = [t]
                pytime = swutil.met2datetime(t, correct=True)
                pystarttime = pytime
            for i in range(len(ulist)):
                if not terse:
                    print("Time(MET + UTCF -> UT):   ", file=pointprint, end="")
                print(
                    "%.3f + %.3f -> %s"
                    % (
                        tlist[i] + 5e-4,
                        ulist[i],
                        swutil.met2fitsString(
                            tlist[i], milliseconds=True, correct=True, format=format
                        ),
                    ),
                    file=pointprint,
                )
                if not terse:
                    print(
                        swutil.met2fitsString(
                            tlist[i],
                            milliseconds=True,
                            correct=True,
                            format="YYMMDD_SOD:                   %y%m%d_%q       DOY=%j",
                        ),
                        file=pointprint,
                    )
            if manualAttitudes:
                pointings = manualAttitudes
            else:
                pointings = PointingEntry.listfortimes([pytime])
            if orbits:
                if debug:
                    debug.write(
                        "About to get the satellite position for %s\n" % (pystarttime,)
                    )
                    debug.flush()
                sat, satpos = orbits.getSatellite(pystarttime)
                satra, satdec, satradius = satpos.radec()
                satlat, satlon = sf_api.wgs84.latlon_of(satpos)
                if debug:
                    debug.write("Got the satellite position\n")
                    debug.flush()
                if terse:
                    print(
                        "Satellite zenith RA=%.2f, dec=%.2f, lat=%.2f N, lon=%.2f E"
                        % (
                            satra._degrees,
                            satdec.degrees,
                            satlat.degrees,
                            satlon.degrees,
                        )
                    )
                else:
                    print(
                        "Swift Zenith(RA,dec):         %.2f, %.2f"
                        % (satra._degrees, satdec.degrees)
                    )
                    print(
                        "Swift Location(lon,lat,alt):  %.2f E, %.2f N, %.0f km"
                        % (
                            satlon.degrees,
                            satlat.degrees,
                            (satradius.km - u.earthRad.to("km")),
                        )
                    )
            for p in pointings:
                if not excelvis:
                    if terse:
                        print(p, file=pointprint)
                    else:
                        print(
                            "Obs Sequence Number:          %08d%03d"
                            % (p._obs, p._segment),
                            file=pointprint,
                        )
                        print(
                            "Obs Target Name:              %s%s"
                            % (p._sourcename, (" (planned)" if p._planned else "")),
                            file=pointprint,
                        )
                        print(
                            "Obs Date and Times:           %s - %s"
                            % (
                                (p._tstart + udelta).strftime("%Y-%m-%d   %H:%M:%S"),
                                (p._tend + udelta).strftime("%H:%M:%S"),
                            ),
                            file=pointprint,
                        )
                        print(
                            "Obs Pointing(ra,dec,roll):    %.3f, %.3f, %.2f"
                            % (p._ra, p._dec, p._roll),
                            file=pointprint,
                        )
                for s in sources:
                    try:
                        s.compute(pystarttime)
                        theta, phi = s.thetangle_phi(p._ra, p._dec, p._roll)
                        # print("exp,cosangle = ",s.exposure(p._ra, p._dec, p._roll) )
                        exp, cosangle = s.exposure(p._ra, p._dec, p._roll)
                        visible = visible or (exp > 0)
                        if orbits:
                            (elev, _, description) = lineofsight(satpos, s)
                            relhoriz = " (%s)" % (description,)
                        else:
                            relhoriz = ""
                        if terse:
                            print(
                                "%s: (theta,phi) = (%.4f, %.4f); exposure = %.0f cm^2%s"
                                % (
                                    s.name,
                                    theta.deg,
                                    phi.wrap_at(180 * u.deg).deg,
                                    exp * cosangle,
                                    relhoriz,
                                ),
                                file=pointprint,
                            )
                        elif excelvis:
                            (elevstart, _, descriptionstart) = lineofsight(
                                orbits.getSatellite(p._tslewend)[1], s
                            )
                            (elevend, _, descriptionend) = lineofsight(
                                orbits.getSatellite(p._tend)[1], s
                            )
                            print(
                                "vis,%s,%s,%s,%s,%.0f,%.1f,%.1f"
                                % (
                                    s.name,
                                    p._tstart.strftime("%Y-%m-%d %H:%M:%S"),
                                    p._tslewend.strftime("%Y-%m-%d %H:%M:%S"),
                                    p._tend.strftime("%Y-%m-%d %H:%M:%S"),
                                    exp * cosangle,
                                    elevstart,
                                    elevend,
                                ),
                                file=pointprint,
                            )
                        else:
                            print(
                                "%s_imageloc (boresight_dist, angle): (%.0f, %.0f)"
                                % (s.name, theta.deg, phi.wrap_at(180 * u.deg).deg),
                                file=pointprint,
                            )
                            print(
                                "%s_exposure (cm^2 cos adjusted): %.0f"
                                % (s.name, exp * cosangle),
                                file=pointprint,
                            )
                            if len(relhoriz) > 3:
                                print(
                                    "%s_altitude: %.0f (%s)"
                                    % (s.name, elev, description),
                                    file=pointprint,
                                )
                    except:
                        traceback.print_exc()
                        # traceback.print_tb(sys.exc_info()[2])
                        print(str(s.catnum) + "ran into trouble in calculations")
                        return
            if visible or not visible_only:
                sys.stdout.write(pointprint.getvalue())
    except getopt.GetoptError:
        usage(argv[0])
    except:
        usage(argv[0])
        traceback.print_exc()
        if debug:
            traceback.print_exc(None, debug)
        # traceback.print_tb(sys.exc_info()[2])
    if debug:
        debug.write("%s: swinfo main() completed\n" % (datetime.datetime.now(),))
        debug.flush()


def checkParse():
    """Go through all the parsed mail messages to see which ones weren't parsed properly"""
    for f in glob.glob("/tmp/latesttrigger.*.processed"):
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
        swinfo_main(sys.argv, debug=debug)
    else:
        swinfo_main([l for l in os.popen("/usr/bin/pbpaste | tr '\r' '\n'")])
