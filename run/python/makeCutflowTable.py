#!/bin/env python

# produce a latex cutflow table from the histograms produced by SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, sys, glob
import ROOT as r

from NavUtils import getAllHistoNames, classifyHistoByName, organizeHistosByType, HistoType, HistoNameClassifier
from SampleUtils import guessSampleFromFilename

r.gROOT.SetBatch(1)

#######
# begin parameters
channel = 'ee'            # should be any of ['ee', 'em', 'mm', 'all']
referenceHisto = 'onebin' # name of the histogram to be used to extract the counts
referenceSyst = 'NOM'     # name of the syst to be used to extract the counts
prodTag = 'Jan15_n0115'   # tag of the files to be used
inputDir = '/export/home/gerbaudo/workarea/Susy2013/SusyTest0/run/anaplots/merged'
# end parameters
#######


inputFileNames = glob.glob(inputDir+'/'+'*'+prodTag+'*.root')
inputFiles = [r.TFile.Open(f) for f in inputFileNames]
assert len(inputFileNames)==len(inputFiles),"Cannot open some of the input files"
print 'input files:\n'+'\n'.join(inputFileNames)

refHistoType = HistoType(pr='', ch=channel, var=referenceHisto, syst=referenceSyst)
histoNameClassifier = HistoNameClassifier()

histosByType = collections.defaultdict(list)

for fname, infile in zip(inputFileNames, inputFiles) :
    samplename = guessSampleFromFilename(fname)
    histoNames = [n for n in getAllHistoNames(infile, onlyTH1=True)
                  if refHistoType.matchAllAvailabeAttrs( histoNameClassifier.histoType( n ) )]
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos : classifyHistoByName(h)
    organizeHistosByType(histosByType, histos, samplename)

refHistos = histosByType # already filtered histonames, all histosByType are refHistos
allSamples = list(set([h.sample for histos in refHistos.values() for h in histos]))
allSelects = sorted(list(set([k.pr for k in histosByType.keys()])))
print 'allSamples : ',allSamples
print 'allSelects : ',allSelects

sampleCountsPerSel = dict() # counts[sample][sel]
countsSampleSel = dict([(s, collections.defaultdict(float)) for s in allSamples])
for t, histos in refHistos.iteritems() :
    if t.ch != channel : continue
    for h in histos :
        sample, sel = h.sample, h.type.pr
        countsSampleSel[sample][sel] += h.Integral()


endrow = ' \\\\'
tablePreamble = '% \usepackage{booktabs}\n' \
                +'% \usepackage{placeins}\n' \
                +'\\begin{table}[htbp] \n' \
                + '\\begin{center} \n' \
                + '\\begin{tabular}{l' + 'r'*len(allSamples) + '} \n' \
                + '\\toprule '
tableEpilogue = '\\bottomrule \n' \
                +'\\end{tabular} \n' \
                +'\\caption{Add a caption here} \n' \
                +'\\end{center} \n' \
                +'\\end{table} \n' \
                +'\\FloatBarrier \n'
colWidth = 12
fwidthField = '%'+str(colWidth)+'s'
header = ' & '.join([fwidthField % t for t in ['selection']+allSamples]) + endrow


print
print tablePreamble
print header
print '\\midrule'
for sel in allSelects : # should really use list comprehension
    counts = [countsSampleSel[sam][sel]
              if sam in countsSampleSel and sel in countsSampleSel[sam] else None \
              for sam in allSamples ]
    line = ' & '.join([fwidthField % f for f in [sel]+["%.1f" % c for c in counts]]) + endrow
    print line
print tableEpilogue

