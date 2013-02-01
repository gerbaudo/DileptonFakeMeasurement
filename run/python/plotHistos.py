#!/bin/env python

# Plot all the histograms produced by SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, sys, glob
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from NavUtils import getAllHistoNames, HistoNameClassifier, organizeHistosByType, setHistoType, setHistoSample
from SampleUtils import colors, guessSampleFromFilename

inputDir = '/export/home/gerbaudo/workarea/Susy2013/SusyTest0/run/anaplots/merged'
inputFileNames = glob.glob(inputDir+'/'+'*.root')
print 'input files:\n'+'\n'.join(inputFileNames)
inputFiles = [r.TFile.Open(f) for f in inputFileNames]


histosByType = collections.defaultdict(list)
classifier = HistoNameClassifier()

for fname, infile in zip(inputFileNames, inputFiles) :
    print '-'*3 + fname + '-'*3
    samplename = guessSampleFromFilename(fname)
    histoNames = getAllHistoNames(inputFiles[0], onlyTH1=True) # [:10] # get only 10 histos for now
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos :
        setHistoType(h, classifier.histoType(h.GetName()))
        setHistoSample(h, samplename)
    organizeHistosByType(histosByType, histos)

def plotHistos(histosDict={'ttbar':None, 'zjets':None},
               extensions=['eps', 'png'], outdir='./plots',
               verbose=False) :
    allHistosEmpty = all([h.GetEntries()==0 for h in histosDict.values()])
    if allHistosEmpty : return
    hnames = [h.GetName() for h in histosDict.values()]
    #if '_pt_' not in hnames[0] : return
    assert 1 == len(set(hnames)),"some histos have different names, something is wrong:\n%s"%str(set(hnames))
    hname = hnames[0]
    if verbose : print "got %d histos for '%s' (samples : %s)" % (len(histosDict), hname, str(histosDict.keys()))
    can = r.TCanvas('can_'+hname, hname, 800, 600)
    can.cd()
    stack = r.THStack('stack_'+hname,'')
    trX, trY = can.GetRightMargin(), can.GetTopMargin()
    legWidth, legHeight = 0.35, 0.35
    leg = r.TLegend(1.0 - trX - legWidth, 1.0 - trY - legHeight, 1.0 - trX, 1.0 - trY)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    firstHisto = None
    for s in ['diboson', 'ttbar', 'zjets', 'multijets'] :
        if s not in histosDict : continue
        h = histosDict[s]
        if not firstHisto : firstHisto = h
        h.SetFillColor(colors[s])
        h.SetDrawOption('bar')
        stack.Add(h)
        leg.AddEntry(h, s+" (%.2f)"%h.Integral(), 'F')
    stack.Draw('hist')
    stack.GetXaxis().SetTitle(firstHisto.GetXaxis().GetTitle())
    stack.GetYaxis().SetTitle(firstHisto.GetYaxis().GetTitle())
    leg.Draw()
    channel, plotRegion = firstHisto.type.ch, firstHisto.type.pr
    def writeLabel(can, label, font='') :
        tex = r.TLatex(0.0, 0.0, '')
        tex.SetNDC()
        if font : tex.SetTextFont(font)
        tex.SetTextAlign(31)
        tex.DrawLatex(1.0-can.GetTopMargin(), 1.0-can.GetRightMargin(), label)
    writeLabel(can, channel+', '+plotRegion, firstHisto.GetTitleFont())
    can.Update()
    for ext in extensions :
        can.SaveAs(outdir+'/'+hname+'.'+ext)

for k,v in histosByType.iteritems() :
    plotHistos(histosDict=dict([(h.sample, h) for h in v]))

               
