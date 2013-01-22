#!/bin/env python

import collections, sys, glob
import ROOT as r

from NavUtils import getAllHistoNames, classifyHistoByName, organizeHistosByType

r.gROOT.SetBatch(1)

inputDir = '/export/home/gerbaudo/workarea/Susy2013/SusyTest0/run/anaplots/merged'
inputFileNames = glob.glob(inputDir+'/'+'*.root')
print 'input files:\n'+'\n'.join(inputFileNames)
inputFiles = [r.TFile.Open(f) for f in inputFileNames]

colors = {
    'ttbar'     : r.kRed+1,
    'zjets'     : r.kOrange-2,
    'wjets'     : r.kBlue-2,
    'diboson'   : r.kSpring+1,
    'singletop' : r.kAzure-4,
    'multijet'  : r.kWhite
    }


def guessSampleFromFilename(filename='', verbose=False) :
    if 'top_' in filename : return 'ttbar'
    elif 'Zjet_' in filename : return 'zjets'
    elif 'ZZ_' in filename \
         or 'WW_' in filename \
         or 'WZ_' in filename : return 'diboson'
    elif 'Wjet_' in filename : return 'wjets'
    else :
        if verbose : print "cannot guess samplename for %s" % filename
    


def exploreFile(file) :
    file.ls()
#exploreFile(inputFiles[0])


histoNames = getAllHistoNames(inputFiles[0], onlyTH1=True)[:10] # get only 10 histos for now
histos = [inputFiles[0].Get(hn) for hn in histoNames]
print "collected histos from %s" % inputFiles[0].GetName()



for h in histos : classifyHistoByName(h)

histosByType = collections.defaultdict(list)


for fname, infile in zip(inputFileNames, inputFiles) :
    samplename = guessSampleFromFilename(fname)
    histoNames = getAllHistoNames(inputFiles[0], onlyTH1=True)[:10] # get only 10 histos for now
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos : classifyHistoByName(h)
    organizeHistosByType(histosByType, histos, samplename)

def plotHistos(histosDict={'ttbar':None, 'zjets':None},
               extensions=['eps', 'png'], outdir='./plots',
               verbose=False) :
    #histosDict = sorted(histosDict.iteritems(), key=lambda (k,v): v.GetEntries()) # need it to stack them?
    hnames = [h.GetName() for h in histosDict.values()]
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
        leg.AddEntry(h, s, 'F')
    stack.Draw('hist')
    stack.GetXaxis().SetTitle(firstHisto.GetXaxis().GetTitle())
    stack.GetYaxis().SetTitle(firstHisto.GetYaxis().GetTitle())
    allHistosEmpty = not stack.GetHistogram().Integral()
    if allHistosEmpty : return
    leg.Draw()
    can.Update()
    for ext in extensions :
        can.SaveAs(outdir+'/'+hname+'.'+ext)

for k,v in histosByType.iteritems() :
    plotHistos(histosDict=dict([(h.sample, h) for h in v]))

               
