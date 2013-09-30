#!/bin/env python

# Generic utility functions for SusyTest0
#
# davide.gerbaudo@gmail.com
# 2013-07-25

import difflib
from functools import wraps
import glob
import os
import re
import subprocess
import unittest

def getCommandOutput(command):
    "lifted from supy (https://github.com/elaird/supy/blob/master/utils/io.py)"
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout,stderr = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}

def filterWithRegexp(stringList, regexp) :
    return [d for d in stringList if re.search(regexp, d)]

def findLatestOneOrTwoRootFiles(dir) :
    files = filter(os.path.isfile, glob.glob(dir + "*.root"))
    files.sort(key=lambda x: os.path.getmtime(x))
    return files[-2:] if len(files)>=2 else files
def findLastRootFile(dir) : return findLatestOneOrTwoRootFiles(dir)[-1]
def guessLatestTagFromLatestRootFiles(dir, debug) :
    "Latest tag if there are at least 2 root files; otherwise empty string"
    files = [f.replace('.root','').rstrip() for f in findLatestOneOrTwoRootFiles(dir)]
    def commonSuffix(ll) : return os.path.commonprefix([l[::-1] for l in ll])[::-1]
    suffix = commonSuffix(files)
    if not len(suffix) : return ''
    def cleanupSampleName(suff) :
        "if the suffix includes part of a sample name, it can be isolated by the '.'"
        lastDot = max(0, suffix.rfind('.'))
        return suffix[suffix.find('_', lastDot):]
    if debug : print "guessLatestTagFromLatestRootFiles: cleanup '%s'"%suffix
    return cleanupSampleName(suffix) if len(suffix) else suffix
def guessMonthDayTag(name) :
    "extract a '_Xxx_yy' tag with a 3-char month and a day"
    match = re.search('(?P<tag>'
                      '_(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)_\d+)',
                      name)
    if match : return match.group('tag')
def guessMonthDayTagFromLastRootFile(dir, debug) :
    lastFile = findLastRootFile(dir)
    if lastFile : return guessMonthDayTag(lastFile)
    elif debug : print "guessMonthDayTagFromLastRootFile: no root files"
def isMonthDayTag(tag) : return guessMonthDayTag(tag)
def commonPrefix(list) : return os.path.commonprefix(list)
def commonSuffix(list) : return os.path.commonprefix([l[::-1] for l in list])[::-1]
def longestCommonSubstring(s1, s2) :
    m = difflib.SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
    return s1[m.a : m.a+m.size]

class Memoize :
    """A class to cache cpu-intensive functions.
    Arguments must be hashable.
    See for example
    http://stackoverflow.com/questions/1988804/what-is-memoization-and-how-can-i-use-it-in-python
    """
    def __init__(self, f) :
        self.f = f
        self.memo = {}
    def __call__(self, *args) :
        if not args in self.memo : self.memo[args] = self.f(*args)
        return self.memo[args]
#
# testing
#
class testGuessMdTag(unittest.TestCase) :
    def testKnownValues(self) :
        knownValues = [('out/foo/ttbar_Aug_23.root',     '_Aug_23'),
                       ('foo/baz_Aug/ttbar_Aug_23.root', '_Aug_23'),
                       ('out/foo/ttbar_Aug_2.root',      '_Aug_2'),
                       ('out/foo/ttbar_August_2.root',      None),
                       ]
        for s, tag in  knownValues :
            gTag = guessMonthDayTag(s)
            self.assertEqual(tag, gTag)

if __name__ == "__main__":
    unittest.main()
