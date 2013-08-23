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
verbose      = options.verbose

assert [susyplot, susysel, fakeprob].count(True)==1,"specify one executable"
scriptDir = 'batchScripts'
template  = ''
template += scriptDir+'/templates/susyPlot.sh.template' if susyplot else ''
template += scriptDir+'/templates/susySel.sh.template'  if susysel else ''
template += scriptDir+'/templates/fakeprob.sh.template' if fakeprob else ''
outdir  = ''
outdir += 'susyplot_out' if susyplot else ''
outdir += 'susysel_out'  if susysel  else ''
outdir += 'fakeprob_out' if fakeprob else ''
if not os.path.isdir(outdir)  : os.mkdir(outdir)
inputTemplate = "filelist/%(sample)s.txt"
outScriptTemplate = scriptDir+'/%(sample)s.sh'
outRootTemplate = "%(outdir)s/%(sample)s_%(tag)s.root"

sampleNames = [d.name for d in datasets if not d.placeholder]
sampleNames = filterWithRegexp(sampleNames, regexp)

def listExists(dset='', flistDir='./filelist') : return os.path.exists(flistDir+'/'+dset+'.txt')
def fillInScriptTemplate(sample, input, output, otherOptions, outScript, scriptTemplate) :
    options  = otherOptions
    options += ' --WH-sample' if 'WH' in sample else ''
    outFile = open(outScript, 'w')
    for line in open(scriptTemplate).readlines() :
        line = line.replace('${inp}', input)
        line = line.replace('${out}', output)
        line = line.replace('${opt}', options)
        line = line.replace('${sample}', sample)
        outFile.write(line)
    outFile.close()

for sample in sampleNames :
    missList, regexUnmatch = not listExists(sample), not re.search(regexp, sample)
    if missList or regexUnmatch:
        msg = "# skipping %s (%s)" % (sample, 'no list' if missList else 'regex unmatch')
        print msg
        continue
    input      = inputTemplate%{'sample':sample}
    output     = outRootTemplate%{'outdir':outdir, 'sample':sample, 'tag':batchTag}
    scriptName = outScriptTemplate%{'sample':sample}
    if overwrite or not os.path.exists(scriptName) :
        fillInScriptTemplate(sample, input, output, otherOptions, scriptName, template)

    cmd = "qsub " \
          "-j oe -V " \
          "-N %(jobname)s " \
          "-o batchlog " \
          " %(scripname)s" \
          % \
          {'jobname':"%s%s"%(sample, batchTag), 'scripname':scriptName}
    print cmd
    if submit :
        out = getCommandOutput(cmd)
        if verbose : print out['stdout']

if not submit : print "This was a dry run; use '--submit' to actually submit the jobs"
