#!/bin/env python

# produce a latex cutflow table from the histograms produced by SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, optparse, sys, glob
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from NavUtils import getAllHistoNames, HistoNameClassifier, organizeHistosByType, HistoType, HistoNameClassifier, setHistoType, setHistoSample
from SampleUtils import guessSampleFromFilename
from PickleUtils import dumpToPickle

#########
# default parameters [begin]
validChannels   = ['ee', 'em', 'mm', 'all']
defaultPickle   = 'counts.pkl'
defaultTag      = 'Jan21_n0115'
defaultHisto    = 'onebin'
defaultRefSyst  = 'NOM'
defaultInputDir = '/export/home/gerbaudo/workarea/Susy2013/SusyTest0/run/anaplots/merged'
# default parameters [end]
#########

parser = optparse.OptionParser()
parser.add_option("-c", "--channel", dest="channel", default=validChannels[0],
                  help="possible channels : %s" % str(validChannels))
parser.add_option("-i", "--input-dir", dest="inputdir", default=defaultInputDir,
                  help="input directory (default '%s')" % defaultInputDir)
parser.add_option("-p", "--pickle", dest="pickle", default=defaultPickle,
                  help="save counts to the specified pikle file (default %s)" % defaultPickle)
parser.add_option("-t", "--tag", dest="tag", default=defaultTag,
                  help="production tag (default '%s')" % defaultTag)
parser.add_option("-H", "--histo", dest="histo", default=defaultHisto,
                  help="histogram to get the counts from (default '%s')" % defaultHisto)
parser.add_option("-S", "--syst", dest="syst", default=defaultRefSyst,
                  help="systematic (default '%s')" % defaultRefSyst)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
channel         = options.channel
inputDir        = options.inputdir
prodTag         = options.tag
referenceHisto  = options.histo
referenceSyst   = options.syst
pickleFile      = options.pickle
verbose         = options.verbose
assert channel in validChannels,"Invalid channel %s (should be one of %s)" % (channel, str(validChannels))
inputFileNames = glob.glob(inputDir+'/'+'*'+prodTag+'*.root')
inputFiles = [r.TFile.Open(f) for f in inputFileNames]
assert len(inputFileNames)==len(inputFiles),"Cannot open some of the input files"

if verbose :
    print "Options:\n" \
          + '\n'.join(["%s : %s" % (o, eval(o))
                       for o in ['channel','inputDir','prodTag', 'referenceHisto', 'referenceSyst']])
    print 'Input files:\n'+'\n'.join(inputFileNames)

# navigate the files and collect the histos
refHistoType = HistoType(pr='', ch=channel, var=referenceHisto, syst=referenceSyst)
histoNameClassifier = HistoNameClassifier()
histosByType = collections.defaultdict(list)
classifier = HistoNameClassifier()

for fname, infile in zip(inputFileNames, inputFiles) :
    samplename = guessSampleFromFilename(fname)
    histoNames = [n for n in getAllHistoNames(infile, onlyTH1=True)
                  if refHistoType.matchAllAvailabeAttrs( histoNameClassifier.histoType( n ) )]
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos :
        setHistoType(h, classifier.histoType(h.GetName()))
        setHistoSample(h, samplename)
    organizeHistosByType(histosByType, histos)
refHistos = histosByType # already filtered histonames, all histosByType are refHistos
allSamples = list(set([h.sample for histos in refHistos.values() for h in histos]))
allSelects = sorted(list(set([k.pr for k in histosByType.keys()])))
print 'allSamples : ',allSamples
print 'allSelects : ',allSelects

# get the counts (adding up what needs to be merged by samplename)
sampleCountsPerSel = dict() # counts[sample][sel]
countsSampleSel = dict([(s, collections.defaultdict(float)) for s in allSamples])
for t, histos in refHistos.iteritems() :
    if t.ch != channel : continue
    for h in histos :
        sample, sel = h.sample, h.type.pr
        countsSampleSel[sample][sel] += h.Integral()

if pickleFile : dumpToPickle(pickleFile, countsSampleSel)

# print the table
endrow = ' \\\\'
colWidth = 12
fwidthField = '%'+str(colWidth)+'s'
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

