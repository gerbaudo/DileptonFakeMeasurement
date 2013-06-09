#!/bin/env python

# Make lists of input root files
#
# Modified version of Matt's script
#
# davide.gerbaudo@gmail.com
# Jan 2013

import optparse, subprocess
import datasets
import re

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

dlist = []
dir = basedir[mode]
cmd = 'ls '+dir+' | grep '+tag+' | grep user'
if verbose : print cmd
ls = subprocess.Popen([cmd], shell=True,stdout=subprocess.PIPE)
dlist  = [l for l in [ll.lstrip().rstrip() for ll in (ls.stdout.read()).split("\n")]
          if l] # skip empty lines

def isDatasetDir(dirname, datasetname):
    "expect the dirname to be smth like *_<samplename>.SusyNt*"
    return re.search('_'+datasetname+'\.SusyNt', dirname)

def makeFile(dataset, name):
    ls = subprocess.Popen(['ls '+dir+dataset+'/*.root* > '+outdir+'/'+name+'.txt'],shell=True)
    ls.wait()

wanted = wantedDsets[mode]
for ds in dlist:
    if verbose : print ds
    for name in wanted:
        if isDatasetDir(ds,name):
            makeFile(ds,name)


