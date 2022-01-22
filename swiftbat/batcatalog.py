"""
BAT Catalog entries as derived from a bcttb.fits file
"""

"""
Copyright (c) 2019, Triad National Security, LLC
All rights reserved.

Copyright 2019. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.

All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1.       Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2.       Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.       Neither the name of Triad National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

from functools import lru_cache
from astropy.io import fits
import numpy as np
from pathlib import Path
from swiftbat import simbadnames, simbadlocation
import re


class BATCatalog():
    filebasename = 'recent_bcttb.fits.gz'
    thisdir = Path(__file__).parent

    def __init__(self, catalogfilename=None):
        self.cattable = self._cattable(catalogfilename=catalogfilename)
        self.makeindices()

    def __getitem__(self, item):
        """
        Index by catalog number, name, or something that resolves to same simplified name
        Returns only the fist match
        :param item:
        :return:
        """
        return self.getall(item)[0]
        # IMPROVEME: Go to SIMBAD and resolve the source, get its position, and find BAT source near that

    def get(self, item, default=KeyError):
        """
        Index by catalog number, name, or something that resolves to same simplified name
        :param item:
        :param default: What to return if nothing matches; reraise the KeyError if that is the default
        :return:
        """
        try:
            return self.__getitem__(item)
        except KeyError:
            if default == KeyError:
                raise
            else:
                return default

    def getall(self, item):
        try:
            return self.bycatnum[int(item)]
        except (ValueError, KeyError):
            pass
        try:
            return self.byname[item]
        except KeyError:
            pass
        try:
            return self.bysimplename[item]
        except KeyError as e:
            try:
                simbadmatches = self.simbadmatch(item)
                if len(simbadmatches):
                    return simbadmatches
            except:
                pass
            raise e  # The KeyError, even though the last failure was by simbadmatch

        # IMPROVEME: Go to SIMBAD and resolve the source, get its position, and find BAT source near that

    def simbadmatch(self, item, tolerance=0.2):
        """
        What rows match the catalogued 
        :param item: 
        :param tolerance: 
        :return: 
        """
        ra, dec = simbadlocation(item)
        return self.positionmatch(ra, dec, tolerance)

    def positionmatch(self, radeg, decdeg, tolerance=0.2):
        rascale = np.cos(np.deg2rad(decdeg))
        tol2 = tolerance ** 2
        raoff = (((self.cattable['RA_OBJ'] - radeg) + 180) % 360 - 180)  # Handle the 360-0 wrap
        dist2 = (raoff * rascale) ** 2 + (self.cattable['DEC_OBJ'] - decdeg) ** 2
        return self.cattable[dist2 < tol2]

    def allnames(self):
        return set([row['NAME'] for row in self.cattable if row['NAME']])

    def makeindices(self):
        self.bycatnum = {}
        self.byname = {}
        self.bysimplename = {}
        for row in self.cattable:
            if not row['CATNUM']:
                continue
            self.bycatnum.setdefault(row['CATNUM'], []).append(row)
            self.byname.setdefault(row['NAME'], []).append(row)
            self.bysimplename.setdefault(self.simplename(row['NAME']), []).append(row)

    @lru_cache(maxsize=0)
    def _cattable(self, catalogfilename=None):
        # Recent catalog file e.g.
        # wget https://heasarc.gsfc.nasa.gov/FTP/swift/data/trend/2021_05/bat/bcatalog/sw03110986012bcttb.fits.gz -O recent_bcttb.fits.gz
        # Currently  2021_09  0654020283
        if catalogfilename is None:
            catalogfilename = self.thisdir.joinpath(self.filebasename)
        catalog = fits.getdata(catalogfilename)
        return catalog

    def simplename(self, name):
        """
        A simple name is the alphanumerics of the name in lower case.
        This folds + and - dec designations, but in practice not a problem
        :param name:
        :return:
        """
        return re.sub('[^a-z0-9]', '', name.lower())
