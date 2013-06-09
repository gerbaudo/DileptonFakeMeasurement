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

import optparse, os
import datasets

validModes = ['mc12', 'susy', 'data',]
validOtherOptions = ['', '--1fb', '--AD']
defaultBatchTag = '_Jan21_n0115'
defaultOtherOptions = ''


parser = optparse.OptionParser()
parser.add_option("-m", "--mode", dest="mode", default=validModes[0],
                  help="possible modes : %s" % str(validModes))
parser.add_option("-t", "--tag", dest="tag", default=defaultBatchTag,
                  help="batch tag (default '%s')" % defaultBatchTag)
parser.add_option("-O", "--other-opt", dest="otherOptions", default=defaultOtherOptions,
                  help="other options for qsub script (default '%s', possible values %s)" \
                  % (defaultOtherOptions, str(validOtherOptions)))
(options, args) = parser.parse_args()
mode = options.mode
batchTag = options.tag
otherOptions = options.otherOptions
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
