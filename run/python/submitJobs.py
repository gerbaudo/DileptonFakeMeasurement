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

import optparse
import os
import re
import datasets
from utils import getCommandOutput

validModes = ['mc12', 'susy', 'data',]
validOtherOptions = ['', '--1fb', '--AD']
defaultBatchTag = '_Jul25_n0145'
defaultOtherOptions = ''


parser = optparse.OptionParser()
parser.add_option("-m", "--mode", dest="mode", default=validModes[0],
                  help="possible modes : %s" % str(validModes))
parser.add_option("-o", "--overwrite", action="store_true", dest="overwrite", default=False,
                  help="overwrite existing batch scripts")
parser.add_option("-O", "--other-opt", dest="otherOptions", default=defaultOtherOptions,
                  help="other options for qsub script (default '%s', possible values %s)" \
                  % (defaultOtherOptions, str(validOtherOptions)))
parser.add_option("-s", "--sample-regexp", dest="samples", default='.*',
                  help="create filelists only for matching samples (default '.*')")
parser.add_option("-S", "--submit", dest="submit", action='store_true', default=False,
                  help="submit jobs (default dry run)")
parser.add_option("-t", "--tag", dest="tag", default=defaultBatchTag,
                  help="batch tag (default '%s')" % defaultBatchTag)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
mode         = options.mode
batchTag     = options.tag
otherOptions = options.otherOptions
overwrite    = options.overwrite
regexp       = options.samples
submit       = options.submit
scriptDir    = 'batchScripts'
verbose      = options.verbose
assert mode in validModes,"Invalid mode %s (should be one of %s)" % (mode, str(validModes))
assert otherOptions in validOtherOptions, "Invalid otherOptions '%s' (should be one of %s)" % (otherOptions, str(validOtherOptions))

datasets = datasets.wantedDsets[mode]


def listExists(dset='', flistDir='./filelist') : return os.path.exists(flistDir+'/'+dset+'.txt')
def fillInScriptTemplate(dataset, suffix, outputfilename,
                         templateFname='batchScripts/batchSub.sh.template') :
    options = ''
    options = options+' --WH-sample' if 'WH' in dataset else options
    outFname= scriptDir+'/'+dataset+'.sh'
    outFile = open(outFname, 'w')
    for line in open(templateFname).readlines() :
        line = line.replace('${inp}', dataset)
        line = line.replace('${out}', dataset+suffix)
        line = line.replace('${opt}', options)
        outFile.write(line)
    outFile.close()

for d in datasets :
    skip = not listExists(d) or not re.search(regexp, d)
    if skip :
        print "# skipping %s" % d
        continue
    scriptName = scriptDir+'/'+d+'.sh'
    if overwrite or not os.path.exists(scriptName) :
        fillInScriptTemplate(d, batchTag, scriptName)

    cmd = "qsub " \
          "-j oe -V " \
          "-v inp=%(d)s,out=%(d)s%(append)s,opt=%(option)s " \
          "-N %(d)s%(append)s " \
          "-o batchlog " \
          " %(sn)s" \
          % \
          {'d':d, 'append':batchTag, 'sn':scriptName, 'option':otherOptions}
    print cmd
    if submit :
        out = getCommandOutput(cmd)
        if verbose : print out['stdout']

if not submit : print "This was a dry run; use '--submit' to actually submit the jobs"
