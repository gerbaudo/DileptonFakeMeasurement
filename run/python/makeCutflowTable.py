#!/bin/env python

# produce a latex cutflow table from the histograms produced by SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, sys, glob
import ROOT as r

from NavUtils import getAllHistoNames, classifyHistoByName, organizeHistosByType
from SampleUtils import guessSampleFromFilename

r.gROOT.SetBatch(1)

referenceHisto = 'onebin' # name of the histogram to be used to extract the counts
referenceSyst = 'NOM'     # name of the syst to be used to extract the counts
inputDir = '/export/home/gerbaudo/workarea/Susy2013/SusyTest0/run/anaplots/merged'
inputFileNames = glob.glob(inputDir+'/'+'*.root')
print 'input files:\n'+'\n'.join(inputFileNames)
inputFiles = [r.TFile.Open(f) for f in inputFileNames]



histosByType = collections.defaultdict(list)

for fname, infile in zip(inputFileNames, inputFiles) :
    samplename = guessSampleFromFilename(fname)
    histoNames = getAllHistoNames(inputFiles[0], onlyTH1=True)
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos : classifyHistoByName(h)
    organizeHistosByType(histosByType, histos, samplename)

refHistos = dict((k, v) for k, v in histosByType.iteritems()
                 if k.var==referenceHisto and k.syst==referenceSyst)
allSamples = list(set([h.sample for histos in refHistos.values() for h in histos]))
allSelects = sorted(list(set([k.pr for k in histosByType.keys()])))

print allSamples
print allSelects

sampleCountsPerSel = dict() # counts[sample][sel]
channel = 'ee' # ['ee', 'em', 'mm', 'all']
countsSampleSel = dict([(s, collections.defaultdict(float)) for s in allSamples])

print "collected all the necessary histograms"
for t, histos in refHistos.iteritems() :
    if t.ch != channel : continue
    for h in histos :
        sample, sel = h.sample, h.type.pr
        countsSampleSel[sample][sel] += h.Integral()

colWidth = 12
endrow = ' \\\\'
fwidthField = '%'+str(colWidth)+'s'
header = ' & '.join([fwidthField % t for t in ['selection']+allSamples]) + endrow
print header

for sel in allSelects :
    counts = [countsSampleSel[sam][sel]
              if sam in countsSampleSel and sel in countsSampleSel[sam] else None
              for sam in allSamples
              ]              
    line = ' & '.join([fwidthField % f for f in [sel]+["%.1f" % c for c in counts]]) + endrow
    print line
