from setuptools import setup

copyright="""
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
 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1.       Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2.       Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.       Neither the name of Triad National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

try:
    with open("README.md", 'r') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = ''

setup(
    name='swiftbat',
    version='0.1.1',
    packages=['swiftbat'],
    package_data={'':['catalog', 'recent_bcttb.fits.gz']},
    url='https://github.com/lanl/swiftbat_python/',
    license='BSD-3-Clause',
    author='David M. Palmer',
    author_email='palmer@lanl.gov',
    description='Routines for dealing with data from BAT on the Neil Gehrels Swift Observatory',
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={'console_scripts': ['swinfo=swiftbat.swinfo:main']},
    install_requires = ['beautifulsoup4', 'pyephem', 'astropy', 'astroquery', 'numpy'],
    classifiers=['Development Status :: 3 - Alpha', 'Intended Audience :: Science/Research', 'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy', ],
    python_requires='>=3.6',
)
