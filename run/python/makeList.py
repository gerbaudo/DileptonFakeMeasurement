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
(options, args) = parser.parse_args()
mode   = options.mode
outdir = options.outdir
regexp = options.samples
tag    = options.tag
assert mode in validModes,"Invalid mode %s (should be one of %s)" % (mode, str(validModes))

# Directory where files are
basedir = {'data' : '/gdata/atlas/ucintprod/SusyNt/data12_'+tag+'/', # data
           'mc12' : '/gdata/atlas/ucintprod/SusyNt/mc12_'+tag+'/',   # mc backgrounds
           'susy' : '/gdata/atlas/ucintprod/SusyNt/susy_'+tag+'/',   # mc signals
           }
tags = [tag]
wantedDsets = datasets.wantedDsets
def filterWithRegexp(stringList, regexp) :
    return [d for d in stringList if re.search(regexp, d)]
wantedDsets = dict([(k, filterWithRegexp(vals, regexp)) for k, vals in wantedDsets.iteritems()])
###############################################################################################
#                           Don't need to edit below here!!!                                  #
###############################################################################################
dlist = []
dir = basedir[mode]
for tag in tags :
    print tag
    print "ls " + dir + " | grep " + tag + " | grep user"
    ls = subprocess.Popen(["ls " + dir + " | grep " + tag + " | grep user"],shell=True,stdout=subprocess.PIPE)
    dlist  = dlist + [l for l in [ll.lstrip().rstrip()
                                  for ll in (ls.stdout.read()).split("\n")]
                      if l] # skip empty lines

def contains(dataset, name):
    if (name + ".") in dataset:
        return True
    return False

def makeFile(dataset, name):
    ls = subprocess.Popen(['ls '+dir+dataset+'/*.root* > '+outdir+'/'+name+'.txt'],shell=True)
    ls.wait()

wanted = wantedDsets[mode]
for ds in dlist:
    print ds
    for name in wanted:
        if(contains(ds,name) and not ("_a" in ds)):
            makeFile(ds,name)


