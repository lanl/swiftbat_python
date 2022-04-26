from __future__ import division, print_function

import os
import glob
import datetime
import re
from pathlib import Path
import ftplib

"""
From the FITS file:
COMMENT
COMMENT  Swift clock correction file
COMMENT
COMMENT  This file is used by tools that correct Swift times for the measured
COMMENT  offset of the spacecraft clock based on fits to data provided by the
COMMENT  Swift MOC. The columns contain polynomial coefficients for various time
COMMENT  ranges and are used to compute the appropriate clock correction for any
COMMENT  given mission elapsed time (after the first post-launch measurement).
COMMENT  The technique is very similar to that used for RXTE fine timing.
COMMENT
COMMENT  In contrast to RXTE where these fine clock corrections have magnitudes
COMMENT  measured in the tens of microseconds, the Swift corrections are much
COMMENT  larger (~0.7 seconds at launch). At this writing (June 2005) the
COMMENT  spacecraft clock still lags behind reference clocks on the ground but
COMMENT  also runs fast.  The spacecraft clock is currently accelerating slowly.
COMMENT
COMMENT  Time is divided into intervals, expressed in spacecraft MET.
COMMENT  Columns are: TSTART TSTOP TOFFSET C0 C1 C2 TYPE CHI_N DOF
COMMENT  where
COMMENT     TSTART and TSTOP are the interval boundaries (MET),
COMMENT     TOFFSET is the clock offset (ie, TIMEZERO),
COMMENT     MAXGAP is the largest gap with no data (-1 if unknown),
COMMENT     CHI_N is the reduced chi-square value (-1 if unknown),
COMMENT     DOF is the number of degrees of freedom,
COMMENT     TYPE indicates the segment clock continuity (0=YES; 1=NO).
COMMENT  (MAXGAP, CHI_N, DOF, and TYPE are retained for continuity with
COMMENT  the file format used by RXTE but the values are currently all
COMMENT  defaults and are not used in computing the actual clock offset.)
COMMENT
COMMENT  For a mission time of T, the correction in seconds is computed
COMMENT  with the following:
COMMENT     T1 = (T-TSTART)/86400
COMMENT     TCORR = TOFFSET + (C0 + C1*T1 + C2*T1*T1)*1E-6
"""

theClockData = None

caveat_time = 86400 * 90  # print a caveat if UTCF is more than 90 days stale


class clockErrData:
    # URL is never used on David Palmer's machine (updates handled by ...swift-trend/getit)
    # clockurl = "https://heasarc.gsfc.nasa.gov/FTP/swift/calib_data/sc/bcf/clock/"
    clockurl = "ftps://heasarc.gsfc.nasa.gov/caldb/data/swift/mis/bcf/clock/"
    clockhost = 'heasarc.gsfc.nasa.gov'
    clockhostdir = '/caldb/data/swift/mis/bcf/clock/'
    clockfile_regex = 'swclockcor20041120v\d*.fits'
    # FIXME this should be derived from the dotswift params
    clocklocalsearchpath = ['/opt/data/Swift/swift-trend/clock',
                            os.path.expanduser('~/.swift/swiftclock'),
                            '/tmp/swiftclock']
    clockfilepattern = 'swclockcor20041120v*.fits'

    def __init__(self):
        try:
            self._clockfile = self.clockfile()
            from astropy.io import fits
            f = fits.open(self._clockfile)
            # copies are needed to prevent pyfits from holding open a file handle for each item
            self._tstart = f[1].data.field('TSTART').copy()
            self._tstop = f[1].data.field('TSTOP').copy()
            self._toffset = f[1].data.field('TOFFSET').copy()
            self._c0 = f[1].data.field('C0').copy()
            self._c1 = f[1].data.field('C1').copy()
            self._c2 = f[1].data.field('C2').copy()
            f.close()
        except:
            self._clockfile = ""
            raise RuntimeError("No clock UTCF file")

    def utcf(self, t, trow=None):  # Returns value in seconds to be added to MET to give correct UTCF
        if not self._clockfile:
            return 0.0, "No UTCF file"
        if trow is None:
            trow = t  # What time to use for picking out the row
        caveats = ""
        # Assume that the times are in order, and take the last one before t
        row = [i for i in range(len(self._tstart)) if self._tstart[i] <= trow]
        if len(row) == 0:
            row = 0
            caveats += "Time before first clock correction table entry"
        else:
            row = row[-1]
            if self._tstop[row] + caveat_time < t:
                caveats += "Time %.1f days after clock correction interval\n" % ((t - self._tstop[row]) / 86400)
        ddays = (t - self._tstart[row]) / 86400.0
        tcorr = self._toffset[row] + 1e-6 * (self._c0[row] + ddays * (self._c1[row] + ddays * self._c2[row]))
        return -tcorr, caveats

    def updateclockfiles(self, clockdir, ifolderthan_days=30, test_days=1):
        """
        Update the clock files if the current clockfile is old and we haven't checked recently for new ones
        :param clockdir:
        :param ifolderthan_days:
        :param test_days:
        :return:
        """
        testfile = os.path.join(clockdir, "clocktest")
        try:
            clockfile = sorted(list(glob.glob(os.path.join(clockdir, self.clockfilepattern))))[-1]
            age = datetime.datetime.utcnow() - datetime.datetime.fromtimestamp(os.path.getmtime(clockfile))
            if age.total_seconds() < (86400 * ifolderthan_days):
                return
            # Check no more than once a day
            testage = datetime.datetime.utcnow() - datetime.datetime.fromtimestamp(os.path.getmtime(testfile))
            if testage.total_seconds() < (86400 * test_days):
                return
        except:
            pass
        # Requires wget.  If this is a problem, use ftplib.FTP
        # os.system(
        #     "wget -q --directory-prefix=%s --no-host --no-clobber --cut-dirs=6 -r %s"
        #     % (clockdir, self.clockurl) )
        try:
            ftps = ftplib.FTP_TLS(self.clockhost)
            ftps.login()    # anonymous
            ftps.prot_p()   # for ftps
            ftps.cwd(self.clockhostdir)
            clockreg = re.compile(self.clockfile_regex)
            ftplatest = sorted([f for f in ftps.nlst() if clockreg.match(f)])[-1]
            locallatest = Path(clockdir).joinpath(ftplatest)
            if not locallatest.exists():
                with open(locallatest, "wb") as newfile:
                    ftps.retrbinary(f'RETR {ftplatest}', newfile.write)
            open(testfile,'w').write(' ') # touch
        except Exception as e:
            print(e, file=sys.stdout)


    def clockfile(self):
        for clockdir in self.clocklocalsearchpath:
            if os.path.exists(clockdir):
                break
        else:
            for clockdir in self.clocklocalsearchpath:
                if os.path.exists(os.path.dirname(clockdir)):
                    # If the parent exists, it is ok to add the clockdir to it if possible,
                    # but don't want to build a whole new directory tree
                    try:
                        os.mkdir(clockdir)
                        break
                    except:
                        pass
            else:
                try:
                    os.makedirs(clockdir)   # Force directory to exist in the temp directory
                except FileExistsError:
                    pass
        self.updateclockfiles(clockdir)
        clockfile = sorted(list(glob.glob(os.path.join(clockdir, self.clockfilepattern))))[-1]
        return clockfile


def utcf(met, printCaveats=True, returnCaveats=False):
    """
    Correction factor to add to MET to get UTC
    :param met:
    :param printCaveats:
    :param returnCaveats:
    :return:
    """
    global theClockData
    if theClockData == None:
        theClockData = clockErrData()
    try:
        uc = [theClockData.utcf(t_) for t_ in met]
        # http://muffinresearch.co.uk/archives/2007/10/16/python-transposing-lists-with-map-and-zip/
        # Not valid after Python 2.7
        u, c = map(None, *uc)
        if printCaveats and any(c):
            print("\n".join(["**** " + c_ for c_ in c if c_]))
    except TypeError:
        u, c = theClockData.utcf(met)
        if printCaveats and c:
            print("**** " + c)
    if returnCaveats:
        return u, c
    else:
        return u
