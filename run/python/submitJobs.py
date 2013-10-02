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
parser.add_option('--fakerate', action='store_true', default=False, help='fake rate with code from Matt')
parser.add_option('--fakepred', action='store_true', default=False, help='run FakePred by Matt')
parser.add_option("-o", "--overwrite", action="store_true", dest="overwrite", default=False,
                  help="overwrite existing batch scripts")
parser.add_option("-O", "--other-opt", dest="otherOptions", default='',
                  help="other options that will be passed on to the executable; double quotes if necessary")
parser.add_option("-s", "--sample-regexp", dest="samples", default='.*',
                  help="create filelists only for matching samples (default '.*')")
parser.add_option("-e", "--exclude-regexp", dest="exclude", default=None,
                  help="submit jobs excluding matching samples (eg. bkg only: '(period|simplified)')")
parser.add_option("-S", "--submit", dest="submit", action='store_true', default=False,
                  help="submit jobs (default dry run)")
parser.add_option("-t", "--tag", dest="tag", default=defaultBatchTag,
                  help="batch tag (default '%s')" % defaultBatchTag)
parser.add_option("--alsoplaceholders", action="store_true", default=False,
                  help="consider dummy samples as well (skip them by default)")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
batchTag     = options.tag
otherOptions = options.otherOptions
overwrite    = options.overwrite
regexp       = options.samples
exclude      = options.exclude
submit       = options.submit
susyplot     = options.susyplot
susysel      = options.susysel
fakepred     = options.fakepred
fakeprob     = options.fakeprob
fakerate     = options.fakerate
alsoph       = options.alsoplaceholders
verbose      = options.verbose

if not [susyplot, susysel, fakeprob, fakerate, fakepred].count(True)==1 :
    parser.error("specify one executable")
scriptDir = 'batchScripts'
template  = ''
template += scriptDir+'/templates/susyPlot.sh.template' if susyplot else ''
template += scriptDir+'/templates/susySel.sh.template'  if susysel else ''
template += scriptDir+'/templates/fakeprob.sh.template' if fakeprob else ''
template += scriptDir+'/templates/fakerate.sh.template' if fakerate else ''
template += scriptDir+'/templates/fakepred.sh.template' if fakepred else ''
outdir = 'out/'
logdir = 'log/'
def subdir() :
    if susyplot : return 'susyplot'
    if susysel  : return 'susysel'
    if fakeprob : return 'fakeprob'
    if fakerate : return 'fakerate'
    if fakepred : return 'fakepred'
def formAndCreateOutdir(basedir, subdir) :
    d = basedir+'/'+subdir
    if not os.path.isdir(d)  : os.makedirs(d)
    return d
outdir = formAndCreateOutdir('out/', subdir())
logdir = formAndCreateOutdir('log/', subdir())
inputTemplate = "filelist/%(sample)s.txt"
outScriptTemplate = scriptDir+'/%(sample)s.sh'
outRootTemplate = "%(outdir)s/%(sample)s_%(tag)s.root"
outLogTemplate = "%(logdir)s/%(sample)s_%(tag)s.log"

sampleNames   = [d.name for d in datasets if not d.placeholder or alsoph]
sampleNames   = filterWithRegexp(sampleNames, regexp)
excludedNames = [] if exclude is None else filterWithRegexp(sampleNames, exclude)
sampleNames   = [s for s in sampleNames if s not in excludedNames]
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
    outlog     = outLogTemplate%{'logdir':logdir, 'sample':sample, 'tag':batchTag}
    scriptName = outScriptTemplate%{'sample':sample}
    fileExists = os.path.exists(scriptName)
    if overwrite or not fileExists :
        fillInScriptTemplate(sample, input, output, otherOptions, scriptName, template)
    elif fileExists : print "warning, not overwriting existing script '%s'"%scriptName
    cmd = "qsub " \
          "-j oe -V " \
          "-N %(jobname)s " \
          "-o %(outlog)s " \
          " %(scripname)s" \
          % \
          {'jobname':"%s%s"%(sample, batchTag), 'outlog':outlog, 'scripname':scriptName}
    print cmd
    if submit :
        out = getCommandOutput(cmd)
        if verbose : print out['stdout']

if not submit : print "This was a dry run; use '--submit' to actually submit the jobs"
