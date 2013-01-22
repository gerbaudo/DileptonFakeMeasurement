#!/bin/env python

import collections, re, sys, glob
import ROOT as r

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

def getAllHistos(inputDir, verbose=False, onlyTH1=False, onlyTH2=False, onlyTH3=False) :
    objectNames = []
    allKeys = [k for k in inputDir.GetListOfKeys()]
    directoryKeys = [k for k in allKeys if r.TClass(k.GetClassName()).InheritsFrom(r.TDirectory.Class())]
    directoryNames = [k.GetName() for k in directoryKeys]
    def isTH1key(key) :
        kc = r.TClass(key.GetClassName())
        return kc.InheritsFrom(r.TH1.Class()) \
               and not kc.InheritsFrom(r.TH2.Class()) \
               and not kc.InheritsFrom(r.TH3.Class())
    def isTH2key(key) : return r.TClass(key.GetClassName()).InheritsFrom(r.TH2.Class())
    def isTH3key(key) : return r.TClass(key.GetClassName()).InheritsFrom(r.TH3.Class())
    if onlyTH1 : allKeys = [k for k in allKeys if isTH1key(k)]
    if onlyTH2 : allKeys = [k for k in allKeys if isTH2key(k)]
    if onlyTH3 : allKeys = [k for k in allKeys if isTH3key(k)]
    histoNames = [k.GetName() for k in allKeys if r.TClass(k.GetClassName()).InheritsFrom(r.TH1.Class())]
    if verbose : print histoNames
    for directory in directoryNames :
        histoNames += getAllHistos(inputDir.Get(directory), verbose, onlyTH1, onlyTH2, onlyTH3)
    return histoNames

histoNames = getAllHistos(inputFiles[0], onlyTH1=True)[:10] # get only 10 histos for now
histos = [inputFiles[0].Get(hn) for hn in histoNames]
print "collected histos from %s" % inputFiles[0].GetName()

class HistoType(object):
    "type of histogram, defined by plot region, channel, variable, syst"
    def __init__(self, pr='', ch='', var='', syst=''):
        for att in ['pr', 'ch', 'var', 'syst'] : setattr(self, att, eval(att))
    def sameas(self, rhs):
        return all([getattr(self,att)==getattr(rhs,att) for att in ['pr', 'ch', 'var', 'syst']])
    def __eq__(self, other) : return self.sameas(other)
    def __str__(self) : return ', '.join(["%s : %s"%(a, getattr(self,a)) for a in ['pr', 'ch', 'var', 'syst']])
    def __hash__(self) : return hash(self.__str__())

def classifyHistoByName(histo, verbose=False) :
    "Extract PR+'_'+chan+'_'+name+'_'+sys from histoname and attach an HistoType attribute"
    n = histo.GetName()
    p = re.compile('(?P<pr>.*?)_'   # plot region (non greedy)
                   +'(?P<ch>.*?)_'  # channel (non greedy)
                   +'(?P<var>.*)_'  # var name (greedy, can contain '_')
                   +'(?P<syst>.*)') # last token, everything that's left
    match = p.search(n)
    if not match :
        if verbose : print "cannot classify %s" % n
        return
    kargs = dict([(g, match.group(g)) for g in ['pr', 'ch', 'var', 'syst']])
    setattr(histo, 'type', HistoType(**kargs))


for h in histos : classifyHistoByName(h)

histosByType = collections.defaultdict(list)

def organizeHistosByType(histosByType = collections.defaultdict(list),
                         histosToOrganize = [], sampleName = '') :
    for h in histosToOrganize :
        setattr(h, 'sample', sampleName)
        histosByType[h.type].append(h)
    return histosByType

for fname, infile in zip(inputFileNames, inputFiles) :
    samplename = guessSampleFromFilename(fname)
    histoNames = getAllHistos(inputFiles[0], onlyTH1=True) #[:10] # get only 10 histos for now
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

               
