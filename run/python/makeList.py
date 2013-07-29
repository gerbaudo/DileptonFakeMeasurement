#!/bin/env python

# Make lists of input root files
#
# davide.gerbaudo@gmail.com
# Jan 2013

import glob
import optparse
import re
import subprocess
from datasets import datasets

validModes = ['mc12', 'susy', 'data',]
defaultTag = 'n0145'

usage="""Create the filelists to be used by SusyPlotter.

Example:
./python/makeList.py -s 'Z(ee|mumu|tautau)' -v
"""
parser = optparse.OptionParser(usage=usage)
parser.add_option("-m", "--mode", dest="mode", default=validModes[0],
                  help="possible modes : %s" % str(validModes))
parser.add_option("-o", "--output-dir", dest="outdir", default='filelist/',
                  help="output directory")
parser.add_option("-s", "--sample-regexp", dest="samples", default='.*',
                  help="create filelists only for matching samples (default '.*')")
parser.add_option("-t", "--tag", dest="tag", default=defaultTag,
                  help="production tag (default '%s')" % defaultTag)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False,
                  help="print even more details, only useful to debug problems")
(options, args) = parser.parse_args()
mode   = options.mode
outdir = options.outdir
regexp = options.samples
tag    = options.tag
verbose= options.verbose
debug  = options.debug
assert mode in validModes,"Invalid mode %s (should be one of %s)" % (mode, str(validModes))
if verbose :
    print "Options:"
    print '\n'.join(["%s : %s" % (o, eval(o))
                     for o in ['mode','outdir','regexp', 'tag',]])

def filterWithRegexp(stringList, regexp) :
    return [d for d in stringList if re.search(regexp, d)]

dsetsNames = [d.name for d in datasets]
dsetsNames = filterWithRegexp(dsetsNames, regexp)

# Directory where files are
basedirs = {'data' : '/gdata/atlas/ucintprod/SusyNt/data12_'+tag+'/', # data
           'mc12' : '/gdata/atlas/ucintprod/SusyNt/mc12_'+tag+'/',   # mc backgrounds
           'susy' : '/gdata/atlas/ucintprod/SusyNt/susy_'+tag+'/',   # mc signals
           }
basedir = basedirs[mode]
dirlist = glob.glob(basedir+'/user.*'+tag)

def isDatasetDir(dirname, datasetname):
    "expect the dirname to be smth like *.<samplename>.SusyNt*"
    return re.search('\.'+datasetname+'\.SusyNt', dirname)
def makeFile(datasetdir, destfilename):
    ls = subprocess.Popen(['ls '+datasetdir+'/*.root* > '+destfilename],shell=True)
    ls.wait()

processedDsets = dict()
nMatching = 0
if verbose : print "Processing..."
for dsdir in dirlist :
    dsname = next((d for d in dsetsNames if isDatasetDir(dsdir, d)), None)
    if not dsname :
        if debug : print "%s dataset not in list...skip it"%dsdir
        continue
    nMatching += 1
    if dsname in processedDsets :
        if debug : print "%s already processed...skip it\n%s"%(dsname,
                                                               '\n'.join([dsdir,
                                                                          processedDsets[dsname]]))
        continue
    if verbose : print dsdir
    flistname = outdir+'/'+dsname+'.txt'
    makeFile(dsdir, flistname)
    processedDsets[dsname] = dsdir
if verbose :
    print "considered %d directories" % len(dirlist)
    print "%d directories matching one of the %d known dataset" % (nMatching, len(dsetsNames))
    print "%s filelists created in %s"%(len(processedDsets), outdir)

