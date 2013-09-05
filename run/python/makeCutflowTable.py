#!/bin/env python

# produce a cutflow table from the histograms produced by SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, optparse, sys, glob
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from CutflowTable import CutflowTable
from NavUtils import getAllHistoNames, HistoNameClassifier, organizeHistosByType, HistoType, HistoNameClassifier, setHistoType, setHistoSample
from SampleUtils import guessGroupFromFilename, isBkgSample
from PickleUtils import dumpToPickle

#########
# default parameters [begin]
validChannels   = ['all', 'ee', 'em', 'mm']
defaultChannel  = validChannels[0]
defaultHisto    = 'onebin'
defaultRefSyst  = 'NOM'
# default parameters [end]
#########

usage="""Using the histograms produced by SusyPlot, print a latex cutflow table.

Examples:
> ./python/makeCutflowTable.py anaplots/merged/  -c all -b -d -s '^sr\d$' # background
> ./python/makeCutflowTable.py anaplots/ -c all  -s '^sr\d$'   # signals
"""
parser = optparse.OptionParser(usage=usage)
parser.add_option("-c", "--channel", default=defaultChannel, help="one of : %s"%str(validChannels))
parser.add_option("-d", "--data", action='store_true', default=False, help="print data")
parser.add_option("-b", "--totbkg", action='store_true', default=False, help="print tot. bkg")
parser.add_option("-r", "--rawcounts", action='store_true', default=False, help="raw rather than weighted")
parser.add_option("-H", "--histo", default=defaultHisto, help="histo from which we count (default '%s')" % defaultHisto)
parser.add_option("-s", "--selregexp", default='.*', help="print only mathing selections (default '.*', any sel); example -s '^sr\d$', see http://www.debuggex.com/r/V82_pzhNDT0ukHMR/1")
parser.add_option("-S", "--syst", default=defaultRefSyst, help="systematic (default '%s')" % defaultRefSyst)
parser.add_option("--csv", default=None, help="save csv to file (default to screen)")
parser.add_option("--tex", default=None, help="save tex to file")
parser.add_option("--pkl", default=None, help="save pickle to file")
parser.add_option("-v", "--verbose", action="store_true", default=False, help="print stuff")
(options, args) = parser.parse_args()
channel         = options.channel
histoname       = options.histo
printData       = options.data
printTotBkg     = options.totbkg
rawcnt          = options.rawcounts
referenceHisto  = options.histo
referenceSyst   = options.syst
selRegexp       = options.selregexp
csvFile         = options.csv
pklFile         = options.pkl
texFile         = options.tex
verbose         = options.verbose

inputs, ext = args, '.root'
if len(inputs) < 1 : parser.error("provide at least one input")
if csvFile and not csvFile.endswith('.csv') : parser.error("csv file must end with 'csv'")
if pklFile and not pklFile.endswith('.pkl') : parser.error("pickle file must end with 'pkl'")
if texFile and not texFile.endswith('.tex') : parser.error("latex file must end with 'tex'")
inputFileNames = [f  for i in inputs for f in glob.glob(i if i.endswith(ext) else i+'/*'+ext)]
assert channel in validChannels,"Invalid channel %s (should be one of %s)" % (channel, str(validChannels))
inputFiles = [r.TFile.Open(f) for f in inputFileNames]

if verbose :
    print "Options:\n" \
          + '\n'.join(["%s : %s" % (o, eval(o))
                       for o in ['channel', 'printData', 'printTotBkg',
                                 'rawcnt', 'referenceHisto', 'referenceSyst','selRegexp',
                                 'pickleFile']])
    print 'Input files:\n'+'\n'.join(inputFileNames)

# navigate the files and collect the histos
referenceType = HistoType(pr='', ch=channel, var=referenceHisto, syst=referenceSyst)
histosByType = collections.defaultdict(list)
classifier = HistoNameClassifier()

histoNames = []
for fname, infile in zip(inputFileNames, inputFiles) :
    sample = guessGroupFromFilename(fname)
    setType, setSample = setHistoType, setHistoSample
    def getType(histoName) : return classifier.histoType(histoName)
    def isRightType(histo) : return referenceType.matchAllAvailabeAttrs(histo.type)
    histonamesCached = len(histoNames)>0
    if not histonamesCached : histoNames = getAllHistoNames(infile, onlyTH1=True, nameStem=histoname)
    histos = filter(isRightType, map(lambda hn :
                                     setSample(setType(infile.Get(hn), getType(hn)), sample),
                                     histoNames))
    if not histonamesCached : histoNames = [h.GetName() for h in histos] # after filtering
    organizeHistosByType(histosByType, histos)
refHistos = histosByType # already filtered histonames, all histosByType are refHistos
allSamples = sorted(list(set([h.sample for histos in refHistos.values() for h in histos])))
allSamples += ['totbkg'] if printTotBkg else []
allSelects = sorted(list(set([k.pr for k in histosByType.keys()])))
if verbose : print 'allSamples : ',allSamples
if verbose : print 'allSelects : ',allSelects

# get the counts (adding up what needs to be merged by samplename)
sampleCountsPerSel = dict() # counts[sample][sel]
countsSampleSel = dict([(s, collections.defaultdict(float)) for s in allSamples])
for t, histos in refHistos.iteritems() :
    if t.ch != channel : continue
    for h in histos :
        sample, sel = h.sample, h.type.pr
        cnt = h.GetEntries() if rawcnt else h.Integral()
        skipIt = not printData and sample=='data'
        countIt = not skipIt
        countsSampleSel[sample][sel] += cnt if countIt else 0.0
        if printTotBkg and isBkgSample(sample) :
            countsSampleSel['totbkg'][sel] += cnt

ct = CutflowTable(allSamples, allSelects, countsSampleSel,
                  isRawCount=rawcnt, selectionRegexp=selRegexp)
csv = ct.csv()
print csv
if csvFile :
    with open(csvFile, 'w') as f : f.write(csv)
if texFile :
    with open(texFile, 'w') as f : f.write(ct.latex())
if pklFile :
    dumpToPickle(pklFile, countsSampleSel)

