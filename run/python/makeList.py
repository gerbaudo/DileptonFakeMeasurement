#!/bin/env python

# Make lists of input root files
#
# Modified version of Matt's script
#
# davide.gerbaudo@gmail.com
# Jan 2013

import datasets
import glob
import optparse
import re
import subprocess

validModes = ['mc12', 'susy', 'data',]
defaultTag = 'n0139'

parser = optparse.OptionParser()
parser.add_option("-m", "--mode", dest="mode", default=validModes[0],
                  help="possible modes : %s" % str(validModes))
parser.add_option("-o", "--output-dir", dest="outdir", default='filelist/',
                  help="output directory")
parser.add_option("-s", "--sample-regexp", dest="samples", default='.*',
                  help="create filelists only for matching samples (default '.*'). Example: Alpgen Z+jets -s '^Z(ee|mumu|tautau).*'")
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

# Directory where files are
basedir = {'data' : '/gdata/atlas/ucintprod/SusyNt/data12_'+tag+'/', # data
           'mc12' : '/gdata/atlas/ucintprod/SusyNt/mc12_'+tag+'/',   # mc backgrounds
           'susy' : '/gdata/atlas/ucintprod/SusyNt/susy_'+tag+'/',   # mc signals
           }
wantedDsets = datasets.wantedDsets
def filterWithRegexp(stringList, regexp) :
    return [d for d in stringList if re.search(regexp, d)]
wantedDsets = dict([(k, filterWithRegexp(vals, regexp)) for k, vals in wantedDsets.iteritems()])

if verbose :
    print "Options:\n" \
          + '\n'.join(["%s : %s" % (o, eval(o)) for o in ['mode','outdir','regexp', 'tag',]])

dir = basedir[mode]
dirlist = glob.glob(dir+'/user.*'+tag)

def isDatasetDir(dirname, datasetname):
    "expect the dirname to be smth like *_<samplename>.SusyNt*"
    return re.search('_'+datasetname+'\.SusyNt', dirname)

def makeFile(datasetdir, destfilename):
    ls = subprocess.Popen(['ls '+datasetdir+'/*.root* > '+destfilename],shell=True)
    ls.wait()

wanted = wantedDsets[mode]
for dsdir in dirlist:
    if verbose : print dsdir
    for dsname in wanted:
        if not isDatasetDir(dsdir, dsname) : continue
        flistname = outdir+'/'+dsname+'.txt'
        makeFile(dsdir, flistname)


