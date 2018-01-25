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
