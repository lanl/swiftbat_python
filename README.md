Swiftbat is a set of Python library routines and command-line utilities that have been developed for the purpose
of retrieving, analyzing, and displaying data from NASA's Swift spacecraft, especially the data from the
Swift Burst Alert Telescope (BAT). Development started before the launch of Swift for private use by the
author and was how he learned Python. As a result, it is not a well-packaged, well-written, coherent library.

All of this data is available from the Swift data archive, but some routines in this library use other access
methods that are not available to the general public. These routines will not be useful except to Swift team members.
This software is provided under the '3-clause BSD license'
(https://opensource.org/licenses/BSD-3-Clause) .
It is provided as-is with no expressed or implied warranty of fitness for any purpose.


If you want to find the exposure of BAT to a point in the FOV, use
```
swiftbat.batExposure(theta, phi)
```
where theta is distance from boresight, phi is angle around the boresight, both in radians.

This package also installs a command-line program 'swinfo' that
tells you Swift Information such as what the MET (onboard-clock)
time is, where Swift was pointing, and whether a specific source
was above the horizon and/or in the field of view.
```
% swinfo 2020-05-05T12:34:56 -o -s "cyg X-1"
Swift Zenith(RA,dec):         232.97, -20.46
Swift Location(lon,lat,alt):  -179.31 E, -20.53 N, 549 km
Time(MET + UTCF -> UT):   610374920.889 + -24.888 -> 2020-05-05T12:34:55.999
YYMMDD_SOD:                   200505_45295.999       DOY=126
Obs Sequence Number:          00033349058
Obs Target Name:              SGR 1935+2154
Obs Date and Times:           2020-05-05   12:32:37 - 13:01:34
Obs Pointing(ra,dec,roll):    293.724, 21.897, 62.77
Cyg_X-1_imageloc (boresight_dist, angle): (14, -133)
Cyg_X-1_exposure (cm^2 cos adjusted): 4230
Cyg_X-1_altitude: 29 (up)
```
Use `swinfo --help` for more details.
