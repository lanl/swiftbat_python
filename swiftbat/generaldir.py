#! /usr/bin/env python

from __future__ import print_function, division, absolute_import

"""
generaldir:
Treat http, ftp, and local files equally as a general directory structure
Requires that the http server be apache with autoindex.
The ftp requirements are that the dir listing begin with the flags and end with the name

David Palmer   palmer@lanl.gov
(C) Copyright 2009 Triad National Security, LLC All rights reserved

"""

import sys
import os
from urllib.request import urlopen
import re
import shutil
import stat


def subTemplate(origstring, keydict,
                allcaps=True):  # if string includes '${FOO}' and keydict['FOO'] defined, substitute
    # The Template class in Python 2.4 would be the better way to do it, except that 2.4 was released yesterday (11/30/04)
    # If allcaps is true, the keys of keydict must be all caps, although the ${foo} need not be
    s = origstring
    found = re.compile(r'''\$\{([^}]+)\}''').findall(s)
    # print(s," found ",found)
    for v in found:
        if allcaps:
            vu = v.upper()
        else:
            vu = v
        try:
            if vu in keydict:
                s = s.replace("${" + v + "}", keydict[vu])
            # print("->",s)
        except:
            print("replace of ${%s} failed" % v)
            raise
    return s


class generalDir:
    _rewildsplitter = re.compile(r'''(?P<unwild>(/*[^$/]+[/]+)*)(?P<firstwild>[^/]*)/*(?P<restwild>.*)''')

    def __init__(self, url):
        self.url = url
        self.lastsubpath = None
        self.lastreadout = None

    def getMatches(self, subpath, reg):
        return re.compile(reg).findall(self.get(subpath))

    def get(self, subpath):
        if self.lastsubpath != subpath or self.lastreadout == None:
            # Cache value so that when we request both files and dirs, only one net trans. made
            self.lastsubpath = subpath
            self.lastreadout = self.openSub(subpath).read().decode('utf-8')
        return self.lastreadout

    def exists(self, subpath=''):
        """ Check if referenced object exists by trying to open it """
        try:
            o = self.openSub(subpath)
            o.close()
            return True
        except:
            return False

    def openSub(self, subpath):
        # print("Opening ",self.url + "/" + subpath)
        try:
            # print("\r"+self.url + "/" + subpath,)
            return urlopen(self.url + "/" + subpath)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print(self.url + "/" + subpath)
            raise

    def makeDirectoryForFile(self, fname):
        try:
            os.makedirs(os.path.dirname(fname))
        except:
            pass  # catch exception thrown when dir already present

    def copyToFile(self, subpath, fname):
        self.makeDirectoryForFile(fname)
        shutil.copyfileobj(self.openSub(subpath), open(fname, "wb"))  # order is src, dest

    def matchPath(self, subpath, wildpath, matchdict=None, regexdict=None):
        """ given a subpath and a wildcard path, with named values and
        regular expressions in matchdict and regexdict, find directories that match
        the pattern and the pattern matched as a list of tuples (path, matchingdict)
        """
        # print("subpath = ",subpath)
        if subpath == None: subpath = ''
        (unwild, firstwild, restwild) = self._splitAndSub(wildpath, matchdict, regexdict)
        # print("join (%s,%s)" % (subpath,unwild))
        subpath = os.path.join(subpath, unwild)  # Add the non-wild stuff to the subpath
        # print("->",subpath)
        if not self.exists(subpath):  # if the unwild stuff doesn't exist then this is a dead end
            # print("Deadend : %s / %s" % (self.url,subpath))
            return []
        if not firstwild:  # if there is no wild part
            assert (not restwild)  # there must not be any more wild part
            # we already checked for existence, so this is a good match
            return [(subpath, matchdict.copy())]
        else:  # There is still some wildness left
            results = []
            matchstring = subTemplate(firstwild, regexdict)
            if -1 != matchstring.find("$"):
                print("Did not substitute variable regular expression in %s" % matchstring)
                raise RuntimeError("Did not substitute variable regular expression in %s" % matchstring)
            dirmatchre = re.compile(matchstring + "/*$")  # optional / allowed for directory
            dirlist = self.dirs(subpath)
            # print("in",subpath,"trying to match string",matchstring,"against",dirlist)
            # print("wildness remaining:",restwild)
            for d in dirlist:
                m = dirmatchre.match(d)
                if m:  # d matches the first wild bit
                    newsubpath = os.path.join(subpath, d)
                    # print("adding",d,"to get",newsubpath)
                    md = matchdict.copy()  # don't disturb original from one recursion level up
                    md.update(m.groupdict())  # add newly discovered variables
                    if restwild:  # more to match
                        results.extend(self.matchPath(newsubpath, restwild, md, regexdict))
                    else:  # nothing more
                        results.append((newsubpath, md))
            if not restwild:  # If there is no more wildness, then match might be a file
                filelist = self.files(subpath)
                filematchre = re.compile(matchstring + "$")  # ending slash is forbidden
                for f in filelist:
                    m = filematchre.match(f)
                    if m:  # file match
                        newsubpath = os.path.join(subpath, f)
                        md = matchdict.copy()  # don't disturb original from one recursion level up
                        md.update(m.groupdict())  # add newly discovered variables
                        results.append((newsubpath, md))
            return results

    def _splitAndSub(self, wildpath, matchdict, regexdict):
        """ returns a (topUNwildpath, firstwildness, restofwildpath) tuple with matches and regexes subbed in """
        wildpath = subTemplate(wildpath, matchdict)  # Do all the matching you can
        if -1 == wildpath.find("$"):
            return (wildpath, '', '')
        m = self._rewildsplitter.match(wildpath)
        # print(wildpath+"-->",m.groupdict())
        if not m or (not m.groupdict()['firstwild'] and m.groupdict()['restwild']):
            # no match if totally non-wild and unslashed.  Not /-terminated, but handled above
            assert (False)
            assert (-1 == wildpath.find("$"))
            return (wildpath, '', '')
        return (m.groupdict()['unwild'], m.groupdict()['firstwild'], m.groupdict()['restwild'])


class httpDir(generalDir):
    """ HTTP specialization for the directory generalization
    Works on apache servers with autoindex generation
    A directory is read by reading from its URL with a trailing /.  (Reading without
    a trailing slash gives a redirection that this lib can't follow.)  In the
    text that comes out, subdirectories and files are represented by links where
    the href is the name (with a trailing slash for subdirs.) with the close quotes,
    close angle brackets, and then a repetition of the name (likewise with the trailing /
    for dirs), except the name is truncated if it is too long.   In Swift, filenames
    can be long, but directory names are short enough to not be truncated.
    
    This can be fixed by adding ?F=0  (fancy listing off) to the end of the URL.  See
    http://httpd.apache.org/docs-2.0/mod/mod_autoindex.html
    """
    # 2008-06-30 added star to " " because some HTTP servers do not put a space in front
    dirmatch = re.compile(r'''(?P<foundname>[^\?"/]+)/"> *(?P=foundname)/''', re.MULTILINE | re.IGNORECASE)
    filematch = re.compile(r'''(?P<foundname>[^\?"/]+)"> *(?P=foundname)''', re.MULTILINE | re.IGNORECASE)
    # And some servers do not evevn understand the /?F=0 non-fancy readout
    filenofancymatch = re.compile(r'''HREF="(?P<foundname>[^\?"/]+)">''', re.MULTILINE | re.IGNORECASE)

    # FIXME allow https
    validurlmatch = re.compile("http://", re.IGNORECASE)

    def __init__(self, url):
        if not self.validurlmatch.search(url):
            raise RuntimeError("Not an http url: %s" % url)
        generalDir.__init__(self, url)
        self.filecache = {}
        self.dircache = {}
        # print("HTTP works on %s", url)

    def files(self, subdir=""):
        try:
            return self.filecache[subdir]
        except:
            files = self.getMatches(subdir + '/?F=0', self.filematch)
            self.filecache[subdir] = files
            return files

    def filesNoFancy(self, subdir=""):  #
        files = self.getMatches(subdir + '/?F=0', self.filenofancymatch)
        return files

    def dirs(self, subdir=""):
        try:
            return self.dircache[subdir]
        except:
            dirs = self.getMatches(subdir + '/?F=0', self.dirmatch)
            self.dircache[subdir] = dirs
            return dirs

    def getFullPath(self, urlrelativepathname):
        return os.path.join(self.url, urlrelativepathname)


class ftpDir(generalDir):
    _filelinematch = re.compile(r'''(?<=^-[r-][w-].[r-][w-].r..\s).+$''', re.MULTILINE | re.IGNORECASE)
    _dirlinematch = re.compile(r'''(?<=^d[r-][w-].[r-][w-].r..\s).+$''', re.MULTILINE | re.IGNORECASE)
    _namematch = re.compile("\S+\s*$")  # last string of nonwhite characters before eol

    def __init__(self, url):
        if not re.compile("ftp://", re.IGNORECASE).search(url):
            raise RuntimeError("Not a valid ftp url: %s" % url)
        generalDir.__init__(self, url)

    def files(self, subdir=""):
        lines = self.getMatches(subdir + '/', self._filelinematch)
        return [self._namematch.findall(l.strip())[0] for l in lines]

    def dirs(self, subdir=""):
        lines = self.getMatches(subdir + '/', self._dirlinematch)
        return [self._namematch.findall(l.strip())[0] for l in lines]


class localDir(generalDir):
    def __init__(self, url):
        # print("Trying %s as local file" % url)
        nofileurl = re.compile(r'''(?<=FILE://)([^>]*)''', re.IGNORECASE).findall(url)
        if len(nofileurl):
            self.dirname = os.path.realpath(nofileurl[0])
        else:
            self.dirname = os.path.realpath(url)
        generalDir.__init__(self, "FILE://" + self.dirname)

    def dirs(self, subdir=""):
        d = os.path.join(self.dirname, subdir)
        # IMPROVEME  use os.scandir for python >= 3.5
        if hasattr(os, 'scandir'):
            # Python >= 3.5
            return [x.name for x in os.scandir('/Volumes/DATA/Swift/swift') if x.is_dir()]
        else:
            allindir = os.listdir(d)
            # FIXME chokes on aliases or links to files that don't exist
            # IMPROVEME  use os.path.isdir
            return [x for x in allindir if stat.S_ISDIR(os.stat(os.path.join(d, x))[stat.ST_MODE])]

    def files(self, subdir=""):
        d = os.path.join(self.dirname, subdir)
        if hasattr(os, 'scandir'):
            # Python >= 3.5
            return [x.name for x in os.scandir('/Volumes/DATA/Swift/swift') if x.is_file()]
        else:
            allindir = os.listdir(d)
            return [x for x in allindir if stat.S_ISREG(os.stat(os.path.join(d, x))[stat.ST_MODE])]

    def openSub(self, subpath):
        open(os.path.join(self.dirname, subpath))

    def exists(self, subpath=''):
        """ Check if referenced object exists by trying to open it """
        try:
            # print("does %s exist?" % os.path.join(self.dirname,subpath))
            mode = os.stat(os.path.join(self.dirname, subpath))[stat.ST_MODE]
        except:
            # print("No")
            # print(os.stat(os.path.join(self.dirname,subpath)))
            return False
        if stat.S_ISREG(mode) or stat.S_ISDIR(mode):
            # print("Yes")
            return True
        else:
            # print("Not as a file or directory")
            return False

    def makedirs(self, subpath):
        """ No analog method or HTTP and FTP genDirectories.
        Make the entire string of directories to the given subpath """
        try:
            os.makedirs(os.path.join(self.dirname, subpath))
        except:
            pass

    def getFullPath(self, urlrelativepathname):
        return os.path.join(self.dirname, urlrelativepathname)


def getDir(url):
    for trialclass in (httpDir, ftpDir, localDir):
        # print("trying ", trialclass)
        try:
            return trialclass(url)
        except:
            pass
    print("Can't open %s" % url)


def dive(url, depth, indent=0):
    d = getDir(url)
    if d:
        dirs = d.dirs()
        files = d.files()
        instring = ("%%%0is" % indent) % " "
        print(instring, "%s:" % url)
        instring += "    "
        if len(files): print(instring, files)
        if depth > 0:
            for subdir in dirs:
                dive(url + "/" + subdir, depth - 1, indent + 4)
        elif len(dirs):
            print(instring, "D:", dirs)


def main():
    regexdict = {'SEQNUM': '(?P<SEQNUM>[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9])',
                 'OBSERVATORY': '(?P<OBSERVATORY>sw)', 'VERSION': '(?P<VERSION>[0-9][0-9][0-9])',
                 'TYPE': '(?P<TYPE>\\W+)', 'CODINGSUFFIXES': '(?P<CODINGSUFFIXES>(.gz|.pgp)*)'}
    for url in ('http://heasarc.gsfc.nasa.gov/FTP/swift/data/',
                'http://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2005_04',
                'http://swift.gsfc.nasa.gov/SDC/data/local/data1/data'):
        print(url)
        archive = getDir(url)
        print("----------------")
        dirpaths = archive.matchPath('', '(sw)?${seqnum}(.${version})?', {}, regexdict)
        print("----------------")
        for (obsdir, d) in dirpaths:
            print(obsdir, d['SEQNUM'], d['VERSION'])
    print("----------------")
    # print(archdir.matchPath("","sw${seqnum}.${version}",{},cache.regexdict))
    # print(archdir.matchPath("","sw${seqnum}.023"+"/"+cache.typepatterns['pob.cat']+'${CODINGSUFFIXES}',{},cache.regexdict))
    # print(archdir.matchPath("",cache.hierdict['DIRPATTERN']+"/"+cache.typepatterns['pob.cat']+'${CODINGSUFFIXES}',{},cache.regexdict))
    if len(sys.argv) <= 1:
        # dive("http://swift.gsfc.nasa.gov/SDC/data/local/data1/data/", 2)
        # dive("file:///tmp", 3)
        # dive("ftp://anonymous:palmer%40lanl.gov@heasarc.gsfc.nasa.gov/swift/data/", 2)
        # print(swiftCache()
        pass
    else:
        dive(sys.argv[1], 2)


if __name__ == "__main__":
    main()
