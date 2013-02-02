#!/bin/env python

# Consider the distribution of one variable, and plot sensitivity Z_n
# as a function of the minimum value of the variable.
#
# Inputs:
# - the root files produced by SusyPlot, for signal and backgrounds
#
# davide.gerbaudo@gmail.com
# Jan 2013


import collections, optparse, os, sys, glob
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from NavUtils import getAllHistoNames, HistoNameClassifier, HistoType, organizeHistosByType, setHistoType, setHistoSample
from SampleUtils import colors, guessSampleFromFilename

#########
# default parameters [begin]
validChannels   = ['all', 'ee', 'em', 'mm',]
defaultTag      = 'Jan21_n0115'
defaultHisto    = 'metrel'
defaultRefSyst  = 'NOM'
defaultSigFile  = './anaplots/wA_noslep_WH_2Lep_2_Jan21_n0115.AnaHists.root'
defaultInputDir = './anaplots/merged'
defaultRegions  = ','.join(["sr%d"%i for i in [6,7,8,9]])
# default parameters [end]
#########

parser = optparse.OptionParser()
parser.add_option("-c", "--channel", dest="channel", default=validChannels[0],
                  help="possible channels : %s" % str(validChannels))
parser.add_option("-i", "--input-dir", dest="inputdir", default=defaultInputDir,
                  help="input directory (default '%s')" % defaultInputDir)
parser.add_option("-t", "--tag", dest="tag", default=defaultTag,
                  help="production tag (default '%s')" % defaultTag)
parser.add_option("-H", "--histo", dest="histo", default=defaultHisto,
                  help="histogram to get the counts from (default '%s')" % defaultHisto)
parser.add_option("-R", "--regions", dest="regions", default=defaultRegions,
                  help="plot regions (default '%s')" % str(defaultRegions))
parser.add_option("-s", "--sig-file", dest="sig", default=defaultSigFile,
                  help="signal file, default : %s" % defaultSigFile)
parser.add_option("-S", "--syst", dest="syst", default=defaultRefSyst,
                  help="systematic (default '%s')" % defaultRefSyst)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
channel         = options.channel
inputDir        = options.inputdir
signalFname     = options.sig
prodTag         = options.tag
referenceHisto  = options.histo
plotRegions     = options.regions.split(',')
referenceSyst   = options.syst
verbose         = options.verbose
assert channel in validChannels,"Invalid channel %s (should be one of %s)" % (channel, str(validChannels))
inputFileNames = glob.glob(inputDir+'/'+'*'+prodTag+'*.root') + glob.glob(signalFname)
inputFiles = [r.TFile.Open(f) for f in inputFileNames]
assert len(inputFileNames)==len(inputFiles),"Cannot open some of the input files"


refHistoType = HistoType(pr='', ch=channel, var=referenceHisto, syst=referenceSyst)
histosByType = collections.defaultdict(list)
classifier = HistoNameClassifier()

for fname, infile in zip(inputFileNames, inputFiles) :
    samplename = guessSampleFromFilename(fname)
    histoNames = [n for n in getAllHistoNames(infile, onlyTH1=True)
                  if refHistoType.matchAllAvailabeAttrs( classifier.histoType( n ) )]
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos :
        setHistoType(h, classifier.histoType(h.GetName()))
        setHistoSample(h, samplename)
    histos = [h for h in histos if h.type.pr in plotRegions]
    organizeHistosByType(histosByType, histos)
refHistos = histosByType # already filtered histonames, all histosByType are refHistos

def isSignal(sampleName) : return 'WH_' in sampleName
allSamples = list(set([h.sample for histos in refHistos.values() for h in histos]))
allBkgNames  = [s for s in allSamples if not isSignal(s)]
sigName = next(s for s in allSamples if isSignal(s))
if verbose : print '\n'.join("%s : %s" % (s,l) for s,l in zip(['bkg','sig'], [str(allBkgNames), sigName]))

bkgHistosByType, sigHistosByType = dict(), dict()
for t,histos in histosByType.iteritems() :
    for h in histos:
        n = h.GetName()
        if isSignal(h.sample) : sigHistosByType[t] = h
        elif t in bkgHistosByType : bkgHistosByType[t].Add(h)
        else : bkgHistosByType[t] = h

def buildHistoSigVsMinThres(bkgHisto, sigHisto, bkgErr=0.2) :
    assert bkgHisto.GetNbinsX()==sigHisto.GetNbinsX(),"need the same binning to build scan"
    h = bkgHisto.Clone(bkgHisto.GetName()+'significance')
    xAx = h.GetXaxis()
    xAx.SetTitle('minimum '+xAx.GetTitle())
    bins = range(1, h.GetNbinsX()+1)[::-1]
    nB, nS = 0., 0.
    for b in bins :
        nB = nB + bkgHisto.GetBinContent(b)
        nS = nS + sigHisto.GetBinContent(b)
        zn = r.RooStats.NumberCountingUtils.BinomialExpZ(nS, nB, bkgErr) if nB else 0.0
        h.SetBinContent(b, zn)
        h.SetBinError(b, 0.)
    return h

for t,hSig in sigHistosByType.iteritems() :
    hBkg = bkgHistosByType[t]
    hZn  = buildHistoSigVsMinThres(hBkg, hSig)
    hZn.SetTitle(str(hBkg.type))
    print "%s : max %.3f Z_n at %.2f" % (str(t), hZn.GetMaximum(), hZn.GetBinCenter(hZn.GetMaximumBin()))
    s = "%s_%s_%s" % (t.pr, t.ch, t.var)
    c = r.TCanvas(s, s, 800, 600)
    c.cd()
    hZn.SetStats(0)
    hZn.Draw()
    pname = s+'.png'
    if os.path.exists(pname) : os.remove(pname)
    c.SaveAs(s+'.png')

    
