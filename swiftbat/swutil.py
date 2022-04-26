#! /usr/bin/env python

from __future__ import print_function, division, absolute_import

"""
swutil
Utilities for dealing with Swift Data
Python datatime objects are UTC, MET is spacecraft MET which needs application of UTCF to get datetime
David Palmer   palmer@lanl.gov

"""

"""
Copyright (c) 2018, Triad National Security, LLC.
All rights reserved.

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

This code was developed using funding from the National Aeronautics and Space Administration (NASA).

"""

import os
import sys

# print __file__
# add the path of this module to the searchpath to let helpers in 
execdir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
if not execdir in map(os.path.abspath, sys.path):
    sys.path.append(execdir)

import re
import datetime
from .clockinfo import utcf
import io

# 2004-12-05T19:43:27
refitstime = re.compile(
    r"(?P<year>[0-9]{4,4})-(?P<month>[0-9]{1,2})-(?P<day>[0-9]{1,2})([T ]+(?P<hour>[0-9]{1,2}):(?P<minute>[0-9]{2,2})(:(?P<second>[0-9]{2,2})([.](?P<fracsecond>[0-9]*))?)?)?",
    re.I)
# 2004:329:12:15:07
redoytime = re.compile(
    r"(?P<year>[0-9]{4,4}):(?P<doy>[0-9]{3,3})(:(?P<hour>[0-9]{2,2}):(?P<minute>[0-9]{2,2})(:(?P<second>[0-9]{2,2})([.](?P<fracsecond>[0-9]*))?)?)?")
# JD2454192.8273
rejdtime = re.compile(r"JD[ ]*24(?P<mytruncjd>[0-9]+(.[0-9]*)?)", re.I)
# MJD14192.5273
remjdtime = re.compile(r"MJD[ ]*(?P<mjd>[0-9]+(.[0-9]*)?)", re.I)
# 140308_13764   IPN seconds of day type format
reIPNsod = re.compile(r"(?P<year>[0-9]{2,4})-?(?P<month>[0-9]{2})-?(?P<day>[0-9]{2})_(T?)(?P<sod>[0-9]{5}(\.[0-9]*)?)")
# General time with multiple options for date and time
# Using slashes to separate ymd is disallowed because it could also
# be m/d/y or d/m/y.  '-' or nothing may be (consistently) used
# ydoy , y-doy or y:doy allowed with y 2 or 4 dig, doy always 3 dig
reday_s = (r"""((?P<year>(\d{2})|(\d{4}))"""  # 2 or 4 digit year
           r"""("""  # either DOY or Month,day
           r"""([:-]?(?P<doy>\d{3}))|"""
           r"""((?P<ymdsep>-?)(?P<month>[0-9]{1,2})(?P=ymdsep)(?P<day>[0-9]{1,2}))))""")
# sod must be 5 integer digits plus optional decimal and fraction
retime_s = (r"[ _]*?(([ _ST](?P<sod>[0-9]{5}(\.[0-9]*)?))"
            r"|([ _T](?P<hour>[0-9]{2})"
            r"(?P<tsep>:?)(?P<minute>[0-9]{2})"
            r"((?P=tsep)(?P<second>[0-9]{2}(\.[0-9]*)?))?))")
reGeneral = re.compile(reday_s + retime_s, re.I)
# 2004:329:12:15:07)

swiftepoch = datetime.datetime(year=2001, month=1, day=1, hour=0, minute=0, second=0)
swiftlaunchtime = datetime.datetime(year=2004, month=11, day=20, hour=17, minute=16)
mjdepoch = datetime.datetime(year=1858, month=11, day=17, hour=0, minute=0, second=0)
# mytruncjd is jd with the 24 stripped off the head, but without the 0.5 day adjustment
mytruncjdepoch = mjdepoch - datetime.timedelta(days=0.5)

fitstimeformat = r"%Y-%m-%dT%H:%M:%S"  # for strftime
yearDOYsecstimeformat = r"%Y_%j_%q"  # %q -> SOD in 00000 (non-standard)


def string2datetime(s, nocomplaint=False, correct=False):
    """
    Convert a string (which may be a Swift MET) to a datetime

    In the case of an MET, is corrected for UTCF if correct==True
    :param s: string with time information
    :param nocomplaint: Don't complain if something is wrong
    :param correct: Include the UTCF when converting an MET to UTC
    :return: UTC
    :rtype: datetime.datetime
    """
    try:
        m = reGeneral.match(s)
        if m:
            mgd = m.groupdict()
            try:
                year = int(mgd['year'])
                if year < 100:
                    # Use 1960-2060 for 2-digit years
                    year += 2000 if (year < 60) else 1900
                if mgd['doy'] is not None:
                    doy = int(mgd['doy'])
                    date = datetime.date(year=year) + datetime.timedelta(days=doy - 1)
                else:
                    month = int(mgd['month'])
                    day = int(mgd['day'])
                    date = datetime.date(year=year, month=month, day=day)
                if mgd['sod'] is not None:
                    addseconds = float(mgd['sod'])
                    tod = datetime.time(0)
                else:
                    hour = int(mgd['hour'])
                    minute = int(mgd['minute'])
                    try:
                        addseconds = float(mgd['second'])
                    except:
                        addseconds = 0.0
                    tod = datetime.time(hour=hour, minute=minute)
                d = (datetime.datetime.combine(date, tod)
                     + datetime.timedelta(seconds=addseconds))
                return d
            except Exception as e:
                print(e)
                print("Time parser failed for {} giving {}".format(s, mgd))
        m = reIPNsod.match(s)
        if m:  # 140308_13764
            mgd = m.groupdict()
            year = int(mgd['year'])
            if year < 100:
                # Use 1960-2060 for 2-digit years
                year += 2000 if (year < 60) else 1900
            month = int(mgd['month'])
            day = int(mgd['day'])
            sod = float(mgd['sod'])
            d = datetime.datetime(year=year, month=month, day=day) + datetime.timedelta(seconds=sod)
            return d
        m = refitstime.match(s)
        if m:  # 2004-12-05T19:43:27
            mgd = m.groupdict()
            if mgd['fracsecond'] == None:
                mgd['microsecond'] = 0
            else:
                mgd['microsecond'] = int(1e6 * float("0." + mgd['fracsecond']))
            for k in mgd.keys():
                if mgd[k]:
                    mgd[k] = int(mgd[k])
                else:
                    mgd[k] = 0
            d = datetime.datetime(year=mgd['year'], month=mgd['month'], day=mgd['day'], hour=mgd['hour'],
                                  minute=mgd['minute'], second=mgd['second'], microsecond=mgd['microsecond'])
            return d
        m = redoytime.match(s)
        if m:  # 2004:329:12:15:07
            mgd = m.groupdict()
            if mgd['fracsecond'] == None:
                mgd['microsecond'] = 0
            else:
                mgd['microsecond'] = int(1e6 * float("0." + mgd['fracsecond']))
            for k in mgd.keys():
                if mgd[k]:
                    mgd[k] = int(mgd[k])
                else:
                    mgd[k] = 0
            d = datetime.datetime(year=mgd['year'], month=1, day=1, hour=mgd['hour'], minute=mgd['minute'],
                                  second=mgd['second'], microsecond=mgd['microsecond'])
            d += datetime.timedelta(days=mgd['doy'] - 1)
            return d
        m = remjdtime.match(s)
        if m:  # MJD 14192.5273
            d = mjdepoch + datetime.timedelta(days=float(m.groupdict()['mjd']))
            return d
        m = rejdtime.match(s)
        if m:  # JD2454192.8273
            d = mytruncjdepoch + datetime.timedelta(days=float(m.groupdict()['mytruncjd']))
            return d
        # None of the patterns match, try treating it as a straight number of seconds
        met = float(s)
        if correct:
            utcf_ = utcf(met, not nocomplaint, False)
            met += utcf_
        return met2datetime(met)
    except:
        if nocomplaint:
            return None
        else:
            print(
                "Invalid time '%s'.  Valid formats: 2004-12-05T19:43:27, 2004:329:12:15:07, 123456789.01234, JD2454192.8273, MJD 14192.5273, 140308_13764" % (
                    s), file=sys.stderr)
            raise


def string2met(s, nocomplaint=False, correct=False):
    d = string2datetime(s, nocomplaint, correct=correct)
    if (d == None):
        return 0
    else:
        return datetime2met(d, correct=correct)


def timedelta2seconds(td):
    return td.days * 86400.0 + (td.seconds * 1.0 + td.microseconds * 1e-6)


def met2datetime(met, correct=False):
    if correct:
        met += utcf(met, False, False)
    return swiftepoch + datetime.timedelta(seconds=met)


def datetime2met(dt, correct=False):
    met = timedelta2seconds(dt - swiftepoch)
    if correct:
        met -= utcf(met, False, False)
    return met


def met2fitsString(met, milliseconds=False, correct=False, format=fitstimeformat):
    d = met2datetime(met, correct=correct)
    if milliseconds:
        ms_string = ".%03d" % min(int(d.microsecond / 1000.0 + 0.5), 999)
    else:
        ms_string = ""
    qspformat = format.split("%q")  # Handle the %q -> seconds of day extension
    if len(qspformat) > 1:
        sod_string = ("%05d" % (((60 * d.hour + d.minute) * 60) + d.second,)) + ms_string
        s = sod_string.join([d.strftime(subformat) for subformat in qspformat])
    else:
        s = d.strftime(format) + ms_string
    return s


def met2mjd(met, correct=False):
    if correct:
        met += utcf(met, False, False)
    deltamjd = swiftepoch - mjdepoch
    mjd = met / 86400.0 + deltamjd.days + deltamjd.seconds / 86400.0
    return mjd


def removeNonAscii(s):
    """ Useful utilitiy to fix up, e.g. email files"""
    # http://stackoverflow.com/questions/1342000/how-to-replace-non-ascii-characters-in-string
    return "".join(i for i in s if ord(i) < 128)


def findInSearchPath(envname, pathlist, basename):
    """
    Look for file: If ${envname} in in environment, return it, else search
    each directory in the path list for a file or directory named basename
    """
    try:
        return os.environ[envname]
    except:
        pass
    if type(pathlist) is str:
        pathlist = pathlist.split(":")
    for p in pathlist:
        try:
            stat = os.stat(os.path.join(p, basename))
            return os.path.join(p, basename)
        except:
            pass
    raise ValueError("Could not find %s in paths : %s" % (basename, pathlist))


class TeeFile(io.TextIOWrapper):
    def __init__(self, *args):
        self._outfiles = args
        self._flush = True

    def write(self, towrite):
        for f in self._outfiles:
            f.write(towrite)
            if self._flush:
                f.flush()


# Decorator for timeout
# http://www.saltycrane.com/blog/2010/04/using-python-timeout-decorator-uploading-s3/
# from http://code.activestate.com/recipes/307871-timing-out-function/
import signal


class TimeoutError(Exception):
    def __init__(self, value="Timed Out"):
        self.value = value

    def __str__(self):
        return repr(self.value)


def timeout(seconds_before_timeout):
    def decorate(f):
        def handler(signum, frame):
            raise TimeoutError()

        def new_f(*args, **kwargs):
            old = signal.signal(signal.SIGALRM, handler)
            signal.alarm(seconds_before_timeout)
            try:
                result = f(*args, **kwargs)
            finally:
                signal.signal(signal.SIGALRM, old)
            signal.alarm(0)
            return result

        new_f.func_name = f.func_name
        return new_f

    return decorate


# Decorator to tee stdout to an open file

def TeeStdoutDecorator(fn, __teefile):
    def inner(*args, **kwargs):
        ostdout = sys.stdout
        try:
            sys.stdout = TeeFile(__teefile, sys.stdout)
            ret = fn(*args, **kwargs)
            sys.stdout = ostdout
        except:
            sys.stdout = ostdout
            raise
        return ret

    return inner


def testTeeStdoutDecorator():
    def testTeeStdoutUndecorated(x):
        print(x)
        print(1 / x)

    testTeeStdout = TeeStdoutDecorator(testTeeStdoutUndecorated, open('/tmp/teetest', 'a'))
    testTeeStdout(5)
    testTeeStdout(0)


# http://stackoverflow.com/questions/132058/showing-the-stack-trace-from-a-running-python-application
# 
def dumponsignal(fname='/tmp/python_trace'):
    import threading, sys, traceback
    def dumpstacks(signal, frame):
        id2name = dict([(th.ident, th.name) for th in threading.enumerate()])
        code = []
        for threadId, stack in sys._current_frames().items():
            code.append("\n# Thread: %s(%d)" % (id2name.get(threadId, ""), threadId))
            for filename, lineno, name, line in traceback.extract_stack(stack):
                code.append('File: "%s", line %d, in %s' % (filename, lineno, name))
                if line:
                    code.append("  %s" % (line.strip()))
        print("\n".join(code))
        open(fname, 'a').write("\n".join(code))

    import signal
    signal.signal(signal.SIGQUIT, dumpstacks)


def main(argv):
    for s in argv:
        print("%-30s -> %23s (corrected) = MET %f" % (
        s, string2datetime(s, False, True), datetime2met(string2datetime(s))))


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        main(['20041225_S12345.67',
              '2004-12-25_12345',
              '20041225_12345',
              '2004-12-25_12345.67',
              '2004-12-05T19:43:27.23',
              '2004-12-25',
              '2004:329:12:15:07.45',
              '123456789.01234',
              '20041225_S12345.67'])
