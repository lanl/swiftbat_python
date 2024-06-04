"""
Skyfield imports, various utilities
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


import skyfield.api as sfapi
import astropy as ap
import numpy as np
import datetime
from functools import lru_cache as __lru_cache
from dateutil.parser import parse as parsedate

load = sfapi.Loader("~/.skyfield", verbose=False)
sfts = load.timescale(builtin=True)


@__lru_cache(0)
def loadsfephem():
    return load("de421.bsp")


@__lru_cache(0)
def __units_in_days(units):
    # prefered units in days
    if units is None:
        return 1.0
    units = units.strip().lower()
    if "days".startswith(units.lower()):
        u = 1.0
    elif "seconds".startswith(units):
        u = 1 / 86400.0
    elif "minutes".startswith(units):
        u = 1 / (24 * 60)
    elif "hours".startswith(units):
        u = 1 / 24
    elif "years".startswith(units):
        u = 365.24217
    else:
        raise NotImplementedError(f"Don't know how to convert {units} to days")
    return u


def __dt_in_days(dt, units=None):
    if np.isscalar(dt):
        if isinstance(dt, (float, int)):
            return dt * __units_in_days(units)
        elif isinstance(dt, datetime.timedelta):
            return dt.total_seconds() / 86400  # Don't do a unit conversion
        elif isinstance(dt, ap.TimeDelta):
            return dt.sec / 86400  # Don't do a unit conversion'
        else:
            raise NotImplementedError(f"Don't know how to convert {type(dt)} to days")
    else:
        return np.array([__dt_in_days(dt_, units) for dt_ in dt])


def _asutc(t: datetime.datetime) -> datetime.datetime:
    try:
        return [_asutc(t_) for t_ in t]
    except:
        pass

    # Convert naive or aware timezone to UTC
    if t.tzinfo is not None and t.tzinfo.utcoffset(None) is not None:
        # timezone-aware.  Convert to UTC
        return t.astimezone(datetime.timezone.utc)
    else:  # Naive time
        # Add tzinfo to a naive time that represents a UTC
        return t.replace(tzinfo=datetime.timezone.utc)


def sftime(
    t, plus=None, units="days", stepby=None, nsteps=None, stepto=None, endpoint=True
):
    """Convert other times to SkyField times

    t can be one or more:
        skyfield.time, astropy.time.Time, float(JD tai), pyephem, or str

    Skyfield times can be arrays.
    Skyfield times can be converted to Python with .utc_datetime()



    Args:
        t (_type_): _description_
        plus (_type_, optional): _description_. Defaults to None.
        units (str, optional): _description_. Defaults to "days".
        stepby (_type_, optional): _description_. Defaults to None.
        nsteps (_type_, optional): _description_. Defaults to None.
        stepto (_type_, optional): _description_. Defaults to None.
    """
    if isinstance(t, str):
        t = [t]
        scalar = True
    else:
        try:
            t[0]
            scalar = False
        except (IndexError, TypeError):
            scalar = True
            t = [t]

    if not isinstance(t[0], sfapi.Time):
        if isinstance(t[0], datetime.datetime):
            t = sfts.from_datetimes(_asutc(t))
        elif isinstance(t[0], ap.time.Time):
            t = sfts.from_astropy(t[0])
        elif isinstance(t[0], float):
            t = sfts.tai_jd(t)
        elif hasattr(t[0], "datetime"):  # Pyephem
            t = sfts.from_datetimes([_asutc(t_.datetime()) for t_ in t])
        elif isinstance(t[0], str):
            t = sfts.from_datetimes([_asutc(parsedate(t_)) for t_ in t])
        else:
            raise NotImplementedError(
                f"Don't know how to change {type(t[0])} into skyfield time"
            )

    t_tai = np.array([t_.tai for t_ in t])

    if plus is not None:
        t_tai += __dt_in_days(plus, units=units)

    nstepkw = (stepby is not None) + (nsteps is not None) + (stepto is not None)
    if nstepkw != 0:
        scalar = False

        if stepby is not None:
            stepby = __dt_in_days(stepby, units=units)
        if stepto is None:
            if len(t_tai) >= 2:
                t_tai, stepto_tai = t_tai[:-1], t_tai[-1]
                nstepkw += 1
        else:
            stepto_tai = sftime(stepto).tai
        if len(t_tai) != 1:
            raise RuntimeError(
                "Only one starting point allowed for stepped time arrays"
            )
        if nstepkw != 2:
            raise RuntimeError("Need two (or zero) of {stepby, stepto, nsteps}")
        if stepby is None:
            dt_tai = np.linspace(
                0, stepto_tai - t_tai[0], num=nsteps, endpoint=endpoint
            )
        elif nsteps is None:
            dt_tai = np.arange(0, stepto_tai - t_tai[0], stepby)
        else:
            dt_tai = np.arange(0, nsteps) * stepby
        t_tai = t_tai[0] + dt_tai

    if scalar and len(t_tai) == 1:
        t_tai = t_tai[0]
    return sfts.tai_jd(t_tai)
