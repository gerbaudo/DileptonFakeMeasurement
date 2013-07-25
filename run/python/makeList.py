#!/bin/env python

# Make lists of input root files
#
# davide.gerbaudo@gmail.com
# Jan 2013

import datasets
import glob
import optparse
import re
import subprocess

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
(options, args) = parser.parse_args()
mode   = options.mode
outdir = options.outdir
regexp = options.samples
tag    = options.tag
verbose= options.verbose
assert mode in validModes,"Invalid mode %s (should be one of %s)" % (mode, str(validModes))
if verbose :
    print "Options:\n"
    print '\n'.join(["%s : %s" % (o, eval(o)) for o in ['mode','outdir','regexp', 'tag',]])

def filterWithRegexp(stringList, regexp) :
    return [d for d in stringList if re.search(regexp, d)]
wantedDsets = datasets.wantedDsets
wantedDsets = dict([(k, filterWithRegexp(vals, regexp)) for k, vals in wantedDsets.iteritems()])

# Directory where files are
basedirs = {'data' : '/gdata/atlas/ucintprod/SusyNt/data12_'+tag+'/', # data
           'mc12' : '/gdata/atlas/ucintprod/SusyNt/mc12_'+tag+'/',   # mc backgrounds
           'susy' : '/gdata/atlas/ucintprod/SusyNt/susy_'+tag+'/',   # mc signals
           }
basedir = basedirs[mode]
dirlist = glob.glob(basedir+'/user.*'+tag)

def isDatasetDir(dirname, datasetname):
    "expect the dirname to be smth like *.<samplename>.SusyNt*"
    return re.search('\.'+datasetname+'\..*SusyNt', dirname)
def makeFile(datasetdir, destfilename):
    ls = subprocess.Popen(['ls '+datasetdir+'/*.root* > '+destfilename],shell=True)
    ls.wait()

datasets = wantedDsets[mode]
processedDsets = []
for dsdir in dirlist :
    dsname = next((d for d in datasets if isDatasetDir(dsdir, d)), None)
    if not dsname : continue # not a directory we're interested in
    assert dsname not in processedDsets, "%s has already been processed "%dsname
    if verbose : print dsdir
    flistname = outdir+'/'+dsname+'.txt'
    makeFile(dsdir, flistname)
    processedDsets.append(dsname)
if verbose : print "%s filelists created in %s"%(len(processedDsets), outdir)

