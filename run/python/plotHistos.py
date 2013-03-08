#!/bin/env python

# Plot all the histograms produced by SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, optparse, sys, glob
#import numpy as np # not available, this hurts.
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

from NavUtils import getAllHistoNames, HistoNameClassifier, organizeHistosByType, setHistoType, setHistoSample
from SampleUtils import colors, guessSampleFromFilename

#########
# default parameters [begin]
defaultTag      = 'Feb21_n0115'
defaultInputDir = './anaplots/merged'
defaultSigFile  = './anaplots/wA_noslep_WH_2Lep_3_Feb21_n0115.AnaHists.root'
defaultSigScale = 10.0
# default parameters [end]
#########

parser = optparse.OptionParser()
parser.add_option("-i", "--input-dir", dest="inputdir", default=defaultInputDir,
                  help="input directory (default '%s')" % defaultInputDir)
parser.add_option("-t", "--tag", dest="tag", default=defaultTag,
                  help="production tag (default '%s')" % defaultTag)
parser.add_option("--test", action="store_true", dest="test", default=False,
                  help="test on a few histograms)")
parser.add_option("-s", "--sig-file", dest="sigFname", default=defaultSigFile,
                  help="signal file (default %s)" % defaultSigFile)
parser.add_option("-S", "--sig-scale", dest="sigScale", default=defaultSigScale,
                  help="signal scale factor (default %.1f)" % defaultSigScale)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
inputDir        = options.inputdir
prodTag         = options.tag
signalFname     = options.sigFname
signalScale     = options.sigScale
justTest        = options.test
verbose         = options.verbose

inputFileNames = glob.glob(inputDir+'/'+'*'+prodTag+'*.root') + glob.glob(signalFname)
print 'input files:\n'+'\n'.join(inputFileNames)
inputFiles = [r.TFile.Open(f) for f in inputFileNames]


histosByType = collections.defaultdict(list)
classifier = HistoNameClassifier()

for fname, infile in zip(inputFileNames, inputFiles) :
    print '-'*3 + fname + '-'*3
    samplename = guessSampleFromFilename(fname)
    histoNames = getAllHistoNames(inputFiles[0], onlyTH1=True)
    histoNames = [h for h in histoNames if any([h.startswith(p) for p in ['sr6', 'sr7', 'sr8', 'sr9']])]
    if justTest : histoNames = histoNames[:10] # just get 10 histos to run quick tests
    histos = [infile.Get(hn) for hn in histoNames]
    for h in histos :
        setHistoType(h, classifier.histoType(h.GetName()))
        setHistoSample(h, samplename)
    organizeHistosByType(histosByType, histos)

def isSignal(sampleName) : return 'WH_' in sampleName

def cumsum(l, leftToRight=True) :
    #return numpy.cumsum(l) # not available ?
    return [sum(l[:i]) for i in range(1,len(l)+1)] if leftToRight \
           else [sum(l[-i:]) for i in range(1,len(l)+1)][::-1]
def mergeOuter(bc, nOuter=2) : # add over/underflow in the first/last bin
    return [sum(bc[:nOuter])] + bc[nOuter:-nOuter] + [sum(bc[-nOuter:])]

def cumSumHisto(histo, leftToRight=True) :
    hCs = histo.Clone(histo.GetName()+'_cs')
    nBinsX = 1+hCs.GetNbinsX() # TH1 starts from 1 (0 underflow, N+1 overflow)
    bc = [hCs.GetBinContent(0)] + [hCs.GetBinContent(i) for i in range(1, nBinsX)] + [hCs.GetBinContent(nBinsX+1)]
    bc = cumsum(mergeOuter(bc), leftToRight)
    tot = bc[-1] if leftToRight else bc[0]
    for i, c in enumerate(bc) :
        hCs.SetBinContent(i+1, c/tot if tot else 0.)
        hCs.SetBinError(i+1, 0.)
    hCs.SetMinimum(0.0)
    hCs.SetMaximum(1.0)
    hCs.SetTitle('')
    hCs.SetFillStyle(0)
    return hCs

def cumEffHisto(histoTemplate, bincontents=[], leftToRight=True) :
    h, bc= histoTemplate, bincontents
    assert h.GetNbinsX()==len(bc),"%d bincontents for %d bins"%(len(bc), h.GetNbinsX())
    h = h.Clone(h.GetName()+'_ce')
    tot = bc[-1] if leftToRight else bc[0]
    for i, c in enumerate(bc) :
        h.SetBinContent(i+1, c/tot if tot else 0.)
        h.SetBinError(i+1, 0.)
    h.SetMinimum(0.0)
    h.SetMaximum(1.0)
    h.SetTitle('')
    h.SetFillStyle(0)
    return h

def cloneAndFillHisto(histo, bincontents=[], suffix='', zeroErr=True) :
    h, bc= histo, bincontents
    assert h.GetNbinsX()==len(bc),"%d bincontents for %d bins"%(len(bc), h.GetNbinsX())
    h = h.Clone(h.GetName()+suffix)
    for i, c in enumerate(bc) :
        h.SetBinContent(i+1, c)
        if zeroErr : h.SetBinError(i+1, 0.)
    return h

def plotCumulativeEfficiencyHisto(pad, h, linecolor=r.kBlack, isPadMaster=True) :
    pad.cd()
    h.SetLineColor(linecolor)
    h.SetLineWidth(2)
    isPadMaster = isPadMaster or 0==len([o for o in pad.GetListOfPrimitives()]) # safety net
    drawOption = 'l' if isPadMaster else 'lsame'
    if isPadMaster :
        xA, yA = h.GetXaxis(), h.GetYaxis()
        xA.SetLabelSize(0)
        xA.SetTitle('')
        yA.SetNdivisions(-201)
        yA.SetTitle('eff')
        yA.SetLabelSize(yA.GetLabelSize()*1.0/pad.GetHNDC())
        yA.SetTitleSize(yA.GetTitleSize()*1.0/pad.GetHNDC())
        yA.SetTitleOffset(yA.GetTitleOffset()*pad.GetHNDC())
        yA.CenterTitle()
    h.Draw(drawOption)
    h.SetStats(0)
    return h
def linearTransform(values, targetRange=[0.0,1.0]) :
    xLoT, xHiT = targetRange[0], targetRange[1]
    xLoO, xHiO = min(values), max(values)
    oriRange, tarRange = (xHiO-xLoO), (xHiT-xLoT)
    return [(xLoT + (x-xLoO)*tarRange/oriRange) if oriRange else 0.0
            for x in values]

def plotZnHisto(pad, h, linecolor=r.kBlack, minY=0.0, maxY=1.0) :
    pad.cd()
    h.SetLineColor(linecolor)
    h.SetLineWidth(2)
    h.SetLineStyle(2)
    bc = [h.GetBinContent(i+1) for i in range(h.GetNbinsX())]
    minZn, maxZn = min(bc), max(bc)
    bc = linearTransform(bc+[0.], [minY, maxY])[:-1] # add one 0 so that the min is at least 0
    for i,b in enumerate(bc) : h.SetBinContent(i+1, b)
    h.Draw('lsame')
    x = h.GetXaxis().GetXmax()
    ax = r.TGaxis(x, minY, x, maxY, minZn, maxZn, 001, "+L")
    ax.SetTitle('Z_{n}')
    ax.CenterTitle()
    ax.SetLabelSize(ax.GetLabelSize()*1.0/pad.GetHNDC())
    ax.SetTitleSize(ax.GetTitleSize()*1.0/pad.GetHNDC())
    ax.SetTitleOffset(ax.GetTitleOffset()*pad.GetHNDC())
    ax.SetLineColor(linecolor)
    ax.SetTitleColor(linecolor)
    ax.SetLabelColor(linecolor)
    ax.Draw()
    return [h, ax]

def binContentsWithUoflow(h) :
    nBinsX = h.GetNbinsX()+1
    return [h.GetBinContent(0)] + \
           [h.GetBinContent(i) for i in range(1, nBinsX)] + \
           [h.GetBinContent(nBinsX+1)]
    
def maxSepVerticalLine(hSig, hBkg, yMin=0.0, yMax=1.0) :
    nxS, nxB = hSig.GetNbinsX(), hBkg.GetNbinsX()
    assert nxS==nxB,"maxSepVerticalLine : histos with differen binning (%d!=%d)"%(nxS,nxB)
    bcS = [hSig.GetBinContent(i) for i in range(1,1+nxS)]
    bcB = [hBkg.GetBinContent(i) for i in range(1,1+nxB)]
    def indexMaxDist(bcS, bcB) :
        return sorted([(i,d) for i,d in enumerate([abs(a-b) for a,b in zip(bcS, bcB)])],
                      key= lambda x : x[1])[-1][0]
    iMax = indexMaxDist(bcS, bcB)
    xPos = hSig.GetBinCenter(iMax+1)
    sep = [abs(a-b) for a,b in zip(bcS, bcB)]
    return r.TLine(xPos, yMin, xPos, yMax)

def plotTopPad(pad, hSig, hBkg) :
    nxS, nxB = hSig.GetNbinsX()+1, hBkg.GetNbinsX()+1  # TH1 starts from 1
    assert nxS==nxB,"maxSepVerticalLine : histos with differen binning (%d!=%d)"%(nxS,nxB)
    pad.SetGridy()
    bcS, bcB = binContentsWithUoflow(hSig), binContentsWithUoflow(hBkg)
    leftToRight = True
    bcLS, bcLB = cumsum(mergeOuter(bcS), leftToRight), cumsum(mergeOuter(bcB), leftToRight),
    leftToRight = False    
    bcRS, bcRB = cumsum(mergeOuter(bcS), leftToRight), cumsum(mergeOuter(bcB), leftToRight)
    zn = r.RooStats.NumberCountingUtils.BinomialExpZ
    bkgUnc = 0.2
    znL = [zn(s, b, bkgUnc) if (b>4.0 and s>0.01) else 0.0 for s,b in zip(bcLS, bcLB)]
    znR = [zn(s, b, bkgUnc) if (b>4.0 and s>0.01) else 0.0 for s,b in zip(bcRS, bcRB)]
    leftToRight = max(znL) >= max(znR)
    zn = znL if leftToRight else znR
    hZn = cloneAndFillHisto(hSig, zn, '_zn')
    hCeS = cumEffHisto(hSig, bcLS if leftToRight else bcRS, leftToRight)
    hCeB = cumEffHisto(hBkg, bcLB if leftToRight else bcRB, leftToRight)
    plotCumulativeEfficiencyHisto(pad, hCeS, r.kRed)
    plotCumulativeEfficiencyHisto(pad, hCeB, r.kBlack, False)
    mark = maxSepVerticalLine(hCeS, hCeB)
    mark.SetLineStyle(2)
    mark.Draw()
    gr =  plotZnHisto(pad, hZn, r.kBlue, 0.0, 1.0) # eff go from 0 to 1
    xAx = hSig.GetXaxis()
    x0, x1 = xAx.GetBinLowEdge(xAx.GetFirst()), xAx.GetBinUpEdge(xAx.GetLast())
    midline = r.TLine(x0, 0.5, x1, 0.5)
    midline.SetLineStyle(2)
    midline.SetLineColor(r.kGray)
    midline.Draw()
    return [hCeS, hCeB, mark, gr, midline]


def plotHistos(histosDict={'ttbar':None, 'zjets':None},
                outdir='./plots', extensions=['png',], # 'eps'],
               verbose=False) :
    allHistosEmpty = all([h.GetEntries()==0 for h in histosDict.values()])
    if allHistosEmpty : return
    hnames = [h.GetName() for h in histosDict.values()]
    #if '_pt_' not in hnames[0] : return
    assert 1 == len(set(hnames)),"some histos have different names, something is wrong:\n%s"%str(set(hnames))
    hname = hnames[0]
    if verbose : print "got %d histos for '%s' (samples : %s)" % (len(histosDict), hname, str(histosDict.keys()))
    can = r.TCanvas('can_'+hname, hname, 800, 600)
    splitFraction = 0.85
    can.cd()
    botPad = r.TPad(can.GetName()+'_bot', 'bot pad', 0.0, 0.0, 1.0, splitFraction, 0, 0, 0)
    botPad.SetTopMargin(0)
    r.SetOwnership(botPad, False)
    can.cd()
    can.Update()
    topPad = r.TPad(can.GetName()+'_top', 'top pad', 0.0, splitFraction, 1.0, 1.0, 0, 0)
    topPad.SetBottomMargin(0)
    r.SetOwnership(topPad, False)

    botPad.Draw()
    botPad.cd()
    stack = r.THStack('stack_'+hname,'')
    #rMarg, lMarg, tMarg = can.GetRightMargin(), can.GetLeftMargin(), can.GetTopMargin()
    rMarg, lMarg, tMarg = botPad.GetRightMargin(), botPad.GetLeftMargin(), botPad.GetTopMargin()
    legWidth, legHeight = 0.325, 0.225
    leg = r.TLegend(1.0 - rMarg - legWidth, 1.0 - tMarg - legHeight, 1.0 - rMarg, 1.0 - tMarg)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    firstHisto = None
    for s in ['diboson', 'ttbar', 'zjets', 'multijets'] :
        if s not in histosDict : continue
        h = histosDict[s]
        if not firstHisto : firstHisto = h
        h.SetFillColor(colors[s])
        h.SetDrawOption('bar')
        stack.Add(h)
        leg.AddEntry(h, s+" (%.2f)"%h.Integral(), 'F')
    if not firstHisto :
        if verbose : print "no bkg histos for %s...continue"%hname
        return
    stack.Draw('hist')
    stack.GetXaxis().SetTitle(firstHisto.GetXaxis().GetTitle())
    stack.GetYaxis().SetTitle(firstHisto.GetYaxis().GetTitle())
    hTot = stack.GetHistogram()
    legOnLeftSide = hTot.GetMaximumBin() > 0.5*hTot.GetNbinsX()
    if legOnLeftSide :
        leg.SetX1(lMarg)
        leg.SetX2(lMarg+legWidth)

    signal = next((h for s,h in histosDict.iteritems() if isSignal(s)), None)
    if signal :
        signal.SetLineColor(r.kRed)
        signal.SetLineWidth(2*signal.GetLineWidth())
        signal.Scale(signalScale)
        signal.Draw('same')
        leg.AddEntry(signal, signal.sample+" (x%.1f, %.2f)"%(signalScale, signal.Integral()), 'L')
        if signal.GetMaximum() > stack.GetMaximum() : stack.SetMaximum(1.1*signal.GetMaximum())
    leg.Draw()
    channel, plotRegion = firstHisto.type.ch, firstHisto.type.pr
    def writeLabel(can, label, font='') :
        tex = r.TLatex(0.0, 0.0, '')
        tex.SetNDC()
        if font : tex.SetTextFont(font)
        tex.SetTextAlign(31)
        tex.DrawLatex(1.0-can.GetTopMargin(), 1.0-can.GetRightMargin(), label)
    label = channel+', '+plotRegion
    leg.SetHeader(label)
    botPad.Update()
    can.Update()
    #topMargin = can.GetTopMargin()
    #topPad = r.TPad(can.GetName()+'_top', 'top pad', 0.0, topMargin, 0.0, 1.0, 0, 0, 0)
    can.cd()
    topPad.Draw()
    topPad.cd()
    hSig, hBkg = signal, stack.GetStack().Last() # 'Last' gives a hist w/ the sum
    graphs = None
    if signal.GetNbinsX() > 2 :
        graphs = plotTopPad(topPad, hSig, hBkg)

#    cS = plotCumulativeEfficiency(topPad, signal, r.kRed)
#    cB = plotCumulativeEfficiency(topPad, stack.GetStack().Last())
#    mark = maxSepVerticalLine(cS, cB)
#    mark.Draw()
    
    topPad.Update()
    can.Update()
    for ext in extensions : can.SaveAs(outdir+'/'+hname+'_lin'+'.'+ext)
#    stack.SetMaximum(5.*stack.GetMaximum())
#    stack.SetMinimum(0.25)
#    can.SetLogy()
#    for ext in extensions : can.SaveAs(outdir+'/'+hname+'_log'+'.'+ext)

for k,v in histosByType.iteritems() :
    plotHistos(histosDict=dict([(h.sample, h) for h in v]))

