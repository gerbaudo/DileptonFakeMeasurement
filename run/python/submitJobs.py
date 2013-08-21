#!/bin/env python

# Submit the SusySelection/SusyPlotter jobs
#
# Converted from Matt's runSP.sh
# Should be run from the directory 'run'
#
# Example usage:
# $ python/submitJobs.py --susyplot -s periodA.physics_Egamma
#
# davide.gerbaudo@gmail.com
# Jan 2013

import optparse
import os
import re
import datasets
from utils import getCommandOutput, filterWithRegexp
from datasets import datasets

defaultBatchTag = '_Jul25_n0145'

parser = optparse.OptionParser()
parser.add_option('--susyplot', action='store_true', default=False)
parser.add_option('--susysel',  action='store_true', default=False)
parser.add_option('--fakeprob', action='store_true', default=False)
parser.add_option("-o", "--overwrite", action="store_true", dest="overwrite", default=False,
                  help="overwrite existing batch scripts")
parser.add_option("-O", "--other-opt", dest="otherOptions", default='',
                  help="other options that will be passed on to the executable; double quotes if necessary")
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
susyplot     = options.susyplot
susysel      = options.susysel
fakeprob     = options.fakeprob
scriptDir    = 'batchScripts'
verbose      = options.verbose

assert [susyplot, susysel, fakeprob].count(True)==1,"specify one executable"
template  = ''
template += 'batchScripts/susyPlot.sh.template' if susyplot else ''
template += 'batchScripts/susySel.sh.template'  if susysel else ''
template += 'batchScripts/fakeprob.sh.template' if fakeprob else ''

dsetsNames = [d.name for d in datasets if not d.placeholder]
dsetsNames = filterWithRegexp(dsetsNames, regexp)

def listExists(dset='', flistDir='./filelist') : return os.path.exists(flistDir+'/'+dset+'.txt')
def fillInScriptTemplate(dataset, suffix, outputfilename, templatefilename) :
    options  = otherOptions
    options += ' --WH-sample' if 'WH' in dataset else ''
    outFname= scriptDir+'/'+dataset+'.sh'
    outFile = open(outFname, 'w')
    for line in open(templatefilename).readlines() :
        line = line.replace('${inp}', "filelist/%s.txt"%dataset)
        line = line.replace('${sample}', dataset)
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
        fillInScriptTemplate(d, batchTag, scriptName, template)

    cmd = "qsub " \
          "-j oe -V " \
          "-N %(jobname)s " \
          "-o batchlog " \
          " %(scripname)s" \
          % \
          {'jobname':"%s%s"%(d,batchTag), 'scripname':scriptName}
    print cmd
    if submit :
        out = getCommandOutput(cmd)
        if verbose : print out['stdout']

if not submit : print "This was a dry run; use '--submit' to actually submit the jobs"
