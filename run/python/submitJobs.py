#!/bin/env python

# Submit the SusySelection/SusyPlotter jobs
#
# Converted from Matt's runSP.sh
# Should be run from the directory 'run'
#
# Example usage:
# $ python/submitJobs.py -m data -s periodA.physics_Egamma
#
# davide.gerbaudo@gmail.com
# Jan 2013

import optparse
import os
import re
import datasets
from utils import getCommandOutput, filterWithRegexp
from datasets import datasets

validOtherOptions = ['', '--1fb', '--AD']
defaultBatchTag = '_Jul25_n0145'
defaultOtherOptions = ''

parser = optparse.OptionParser()
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
batchTag     = options.tag
otherOptions = options.otherOptions
overwrite    = options.overwrite
regexp       = options.samples
submit       = options.submit
scriptDir    = 'batchScripts'
verbose      = options.verbose
assert otherOptions in validOtherOptions, "Invalid otherOptions '%s' (should be one of %s)" % (otherOptions, str(validOtherOptions))

dsetsNames = [d.name for d in datasets if not d.placeholder]
dsetsNames = filterWithRegexp(dsetsNames, regexp)

def listExists(dset='', flistDir='./filelist') : return os.path.exists(flistDir+'/'+dset+'.txt')
def fillInScriptTemplate(dataset, suffix, outputfilename,
                         templateFname='batchScripts/susyPlot.sh.template') :
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

for d in dsetsNames :
    missList, regexUnmatch = not listExists(d), not re.search(regexp, d)
    if missList or regexUnmatch:
        print "# skipping %s (%s)" % (d, 'no list' if missList
                                      else ('regex unmatch' if regexUnmatch
                                            else ''))
        continue
    scriptName = scriptDir+'/'+d+'.sh'
    if overwrite or not os.path.exists(scriptName) :
        fillInScriptTemplate(d, batchTag, scriptName)

    cmd = "qsub " \
          "-j oe -V " \
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
