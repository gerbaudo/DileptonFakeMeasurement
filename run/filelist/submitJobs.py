#!/bin/env python

# Submit the SusySelection/SusyPlotter jobs
#
# Converted from Matt's runSP.sh
# Should be run from the directory 'run'
# (todo: move to python and adjust path names)
#
# Example usage:
# $ ./filelist/submitJobs.py > submit.sh
# $ source submit.sh
#
# davide.gerbaudo@gmail.com
# Jan 2013

import os
import datasets

validModes = ['data', 'mc12', 'susy']
validOtherOptions = ['', '--1fb', '--AD']

mode = 'susy'
otherOptions = ''
batchTag = '_Jan21_n0115'

assert mode in validModes,"Invalid mode %s (should be one of %s)" % (mode, str(validModes))
assert otherOptions in validOtherOptions, "Invalid otherOptions '%s' (should be one of %s)" % (otherOptions, str(validOtherOptions))

datasets = datasets.wantedDsets[mode]

def listExists(dset='', flistDir='./filelist') : return os.path.exists(flistDir+'/'+dset+'.txt')

for d in datasets :
    if not listExists(d) :
        print "# skipping %s" % d
        continue
    cmd = "qsub " \
          "-j oe -V " \
          "-v inp=%(d)s,out=%(d)s%(append)s,opt=%(option)s " \
          "-N %(d)s_%(append)s " \
          "-o batchlog " \
          "batchSPSub.sh" \
          % \
          {'d':d, 'append':batchTag, 'option':otherOptions}
    print cmd
