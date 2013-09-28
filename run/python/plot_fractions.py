#!/bin/env python

# Script to compare the fake matrix percentages from Matt and from Davide
#
# Details:

# This scripts takes as input two text files containing a printout of
# the fractions (i.e. samples compositions) printed out from
# FinalNewFake::buildElectronRateSR and FinalNewFake::buildMuonRateSR.
# See for example fc323331d for the formatting implementation, which
# here is specified in 'Entry'.
#
# For each category
# [mu,el] x [conv,qcd,real] x [ttbar,wjet,zjet,dib,bbar]
# the values from the two files and their ratio are plotted
# vs. selection region.
#
# davide.gerbaudo@gmail.com
# 2013-09-19

import os
import ROOT as r

r.gROOT.SetBatch(1)
r.gErrorIgnoreLevel=r.kError # disable root warnings 
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

from utils import commonPrefix, commonSuffix
from rootUtils import (unitLineFromFirstHisto,
                       firstHisto,
                       drawLegendWithDictKeys
                       )
#___________________________________________________________  
class Entry :
    """
    information stored in one block such as

    Fake rate: percentages for sr 'SRmT2a', muon
                    ttbar	 Wjet	 Zjet	 dib	 bbbar
    mu_percent_qcd: 0.000606924, 0.0915346, 0.00646265, 0.000635465, 0.90076, 
    mu_percent_conv: 0, 0, 0, 0, 0, 

    """
    def __init__(self, percType='', srName='', lep='') :
        assert percType in ['Fake rate', 'Real eff'],"invalid ptype %s"%percType
        self.percType = percType
        self.srName = srName
        self.lep = lep
        self.valsPerVar = {}
    def addLine(self, keys, values, varname='mu_percent_qcd') :
        assert len(keys)==len(values),"invalid line\n%s\n%s"%(str(keys), str(values))
        assert varname not in self.valsPerVar, "%s already there"%varname
        self.valsPerVar[varname] = dict(zip(keys, values))
#___________________________________________________________
def parseFile(filename='') :
    def isInteresting(l) :
        tags = ['percentages for', 'ttbar', '_percent_',]
        return any(t in l for t in tags)
    def startsBlock(l) :
        """Test+parse a line such as
        Fake rate: percentages for sr SRmT2a, muon
        """
        if not ': percentages for ' in l : return
        l = l.replace(' percentages for ','')
        words = l.replace(':',',').split(',')
        assert len(words)==3,"expected 3 tokens, from '%s' got %d : %s"%(l, len(words), str(words))
        pType, sr, lep = words[0].strip(), words[1].replace("'","").strip(), words[2].strip()
        return pType, sr, lep
    def labelsSamples(l) :
        """Test for a line such as:
        ttbar    Wjet    Zjet    dib     bbbar
        """
        if not 'ttbar' in l : return
        else :
            labels = l.split()
            assert len(labels)==5,"not a labels line '%s'"%l
            return labels
    def hasValues(l) :
        """Test+parse a line such as:
        mu_percent_qcd:  0.000425829 0.0915684 0.00646507 0.000451178 0.901089
        or
        mu_percent_qcd: 0.000606924, 0.0915346, 0.00646265, 0.000635465, 0.90076,
        """
        vnames = ["%s_percent_%s"%(lep,t) for lep in ['el','mu'] for t in ['qcd','conv','real']]
        if not any(v in l for v in vnames) : return
        words = l.split()
        assert len(words)==6,"invalid values line '%s'"%l
        vname = words[0].replace(':','')
        values = [float(w.replace(',','')) for w in words[1:]]
        return vname, values
    # actually do the job
    entries = []
    entry = None
    lastLabels = None
    for l in open(filename).readlines() :
        l = l.strip()
        if not isInteresting(l) : continue
        start, labels, values = startsBlock(l), labelsSamples(l), hasValues(l)
        if start :
            if entry is not None : entries.append(entry)
            pType, sr, lep = start
            entry = Entry(pType, sr, lep)
        elif labels :
            lastLabels = labels
        elif values :
            vname, values = values
            entry.addLine(lastLabels, values, vname)
        else :
            print "line not parsed '%s'"%l
    return entries
#___________________________________________________________
def resetErrors(h) :
    for b in range(1, h.GetNbinsX()+1) : h.SetBinError(b, 0.0)
def sortLabels(h) : h.LabelsOption('av') # alpha-sorted, vertical
def getLabels(h) : return [h.GetXaxis().GetBinLabel(b) for b in range(1, h.GetNbinsX())]    

def buildHistos(entries=[], regions=[], hprefix='') :
    "assumes entries are all from the same lep and percentageType (and same user)"
    if not len(entries) : return dict()
    allPtypes, allLeps = set([e.percType for e in entries]), set([e.lep for e in entries])
    assert len(allPtypes)==1, "multiple ptypes %s"%str(allPtypes)
    assert len(allLeps)==1, "multiple leptons %s"%str(allLeps)
    regions = regions if regions else list(set([e.srName for e in entries]))
    firstEntry = entries[0] # use it to deduct samples and percTypes
    varTypes = firstEntry.valsPerVar.keys()
    samples = sorted(list(set([s for vals in firstEntry.valsPerVar.values() for s in vals.keys()])))
    hp = hprefix
    histosPerVtypePerSample = {}
    for vt in varTypes :
        hp = "h_%s_%s"%(hprefix, vt)
        histosPerVtypePerSample[vt] = dict(zip(samples,
                                               [r.TH1F("%s_%s"%(hp,s), "%s %s"%(hp,s), 1,0.0,1.0)
                                                for s in samples]))            
    for e in entries :
        def srMatt2srDavide(sr) : return 'sr CR_WHSS' if sr=='sr SRDavide' else sr
        sr = srMatt2srDavide(e.srName)
        for vname, values in e.valsPerVar.iteritems() :
            for sample, val in values.iteritems() :
                histosPerVtypePerSample[vname][sample].Fill(sr,val)
    [resetErrors(h) for hpt in histosPerVtypePerSample.values() for h in hpt.values()]
    [sortLabels(h) for hpt in histosPerVtypePerSample.values() for h in hpt.values()]
    return histosPerVtypePerSample

def buildRatioHistos(histosNum={}, histosDen={}) :
    "assume that the histos are organized in two dict[vt][sample] provided by buildHistos"
    histosPerVtypePerSample = {}
    def sameLists(l1=[], l2=[]) : return len(l1)==len(l2) and sorted(l1)==sorted(l2)
    assert sameLists(histosNum.keys(), histosDen.keys()),"num and den w/ different vtype keys"
    result = {}
    for vt in histosDen.keys() :
        resultPerSample = {}
        hvn, hvd = histosNum[vt], histosDen[vt]
        assert sameLists(hvn.keys(), hvd.keys()),"num and den w/ different sample keys"
        for s in hvd.keys() :
            hnum, hden = hvn[s], hvd[s]
            hNames = [h.GetName() for h in hnum, hden]
            pre, suf = commonPrefix(hNames), commonSuffix(hNames)
            hName = pre+'_ratio_'+suf if len(pre) and len(suf) else '_over_'.join(hNames)
            h = hnum.Clone(hName)
            h.Divide(hnum, hden) # root already checks that the histos have the same bin labels
            resetErrors(h)
            resultPerSample[s] = h
        result[vt] = resultPerSample
    return result

def drawUnitLine(pad, histos) :
    unitLine = unitLineFromFirstHisto(histos)
    unitLine.Draw()
    pad._line = unitLine
    pad.Update()
def getPad(can, splitFraction, top=False, bot=False) :
    assert [top,bot].count(True)==1,"either top or bot"
    pad = can.GetPad(1 if top else 2) # already there?
    pad = pad if pad else r.TPad(can.GetName()+('_top' if top else '_bot'),
                                 'top pad' if top else 'bot pad', 
                                 0.0, splitFraction if top else 0.0, 
                                 1.0, 1.0 if top else splitFraction,
                                 0, 0)
    r.SetOwnership(pad, False)
    if top : pad.SetBottomMargin(0)
    else :
        pad.SetTopMargin(0)
        pad.SetBottomMargin(2.5*pad.GetBottomMargin()) # make space for text labels
    return pad
def draw(pad, histos={}, colors={}, markers={}, lineStyle=1, asFirst=True) :
    pad.cd()
    opt = 'p' if asFirst else 'p same'
    for sample in histos.keys() :
        h, c, m = histos[sample], colors[sample], markers[sample]
        h.SetLineColor(c)
        h.SetLineStyle(lineStyle)
        h.SetMarkerColor(c)
        h.SetMarkerStyle(m)
        h.Draw(opt)
        opt = opt if 'same' in opt else opt+' same'
    pad.Update()

def plotHistos(histos1=[], histos2=[], histosRatio=[], label1='label1', label2='label2',
               canvas=None, canPrefix='', scaleXlabelsBot=1.0,
               colors={}, markers={},lineStyle1=1, lineStyle2=2, splitFraction=0.80) :
    "Assume again that the histos are in dict[sample]"
    can = canvas if canvas is not None else r.TCanvas(canPrefix+'percent_histos','')
    can.Clear()    
    can.cd()
    padTop, padBot = getPad(can, splitFraction, top=True), getPad(can, splitFraction, bot=True)
    # top pad with ratios
    can.cd()
    padTop.Draw()
    padTop.cd()
    markersr = dict(zip(colors.keys(), len(colors.keys())*[r.kOpenDiamond]))
    for h in histosRatio.values() :
        yAx = h.GetYaxis()
        yAx.SetRangeUser(0.0, 2.0)
        yAx.SetNdivisions(-202)
        yAx.SetLabelSize(1.0/(1.0-splitFraction)*yAx.GetLabelSize())
        h.SetTitle('')
        h.SetStats(0)
    #print 'fh: ',firstHisto(histos1).GetName(),' nratios :',len(histosRatio) # debug
    def bc(h) : return ["%.2f"%h.GetBinContent(b) for b in xrange(1, h.GetNbinsX()+1)]
    # print '\n'.join(["%s : %s"%(s, bc(h)) for s,h in histosRatio.iteritems()]) # debug 
    draw(padTop, histosRatio, colors, markersr, asFirst=True)
    drawUnitLine(padTop, histosRatio)
    labelR = "%s/%s ratio"%(label1[0], label2[0])
    drawLegendWithDictKeys(padTop, {label1 : firstHisto(histos1),
                                    label2 : firstHisto(histos2),
                                    labelR : firstHisto(histosRatio)
                                    },
                            legWidth=0.225, legHeight=0.75)
    # bot pad with histos
    can.cd()    
    padBot.Draw()
    padBot.cd()
    markers1 = dict(zip(colors.keys(), len(colors.keys())*[r.kOpenTriangleUp]))
    markers2 = dict(zip(colors.keys(), len(colors.keys())*[r.kOpenTriangleDown]))
    for h in histos1.values()+histos2.values() :
        h.GetYaxis().SetRangeUser(0.0, 1.08)
        h.GetXaxis().SetLabelSize(scaleXlabelsBot*h.GetXaxis().GetLabelSize())
        h.SetStats(0)
    draw(padBot, histos1, colors, markers1, lineStyle1)
    draw(padBot, histos2, colors, markers2, lineStyle2, asFirst=False)
    drawLegendWithDictKeys(padBot, histos1, legWidth=0.225)
    padBot.Update()
    return can
    
#___________________________________________________________
        
if __name__=='__main__' :
    ioDir = 'out/fakerate/'
    fnameDavide = ioDir+'percentages_davide.txt'
    fnameMatt   = ioDir+'percentages_matt.txt'
    entriesD = parseFile(fnameDavide)
    entriesM = parseFile(fnameMatt)
    indent = 0*'  '
    print indent+"any : D[%s] M[%d]"%(len(entriesD), len(entriesM))
    colors, markers = None, None
    for l in ['muon', 'electron'] :
        indent = 1*'  '
        edl = [e for e in entriesD if e.lep==l]
        eml = [e for e in entriesM if e.lep==l]
        print indent+"%s : D[%d] M[%d]"%(l, len(edl), len(eml))
        for pType in ['Fake rate', 'Real eff'] :
            indent = 2*'  '
            print indent+"pType '%s'"%pType
            edlp = [e for e in edl if e.percType==pType]
            emlp = [e for e in eml if e.percType==pType]
            print indent+"%s : D[%d] M[%d]"%(pType, len(edlp), len(emlp))
            regions = list(set([e.srName for e in edlp+emlp]))
            samples = sorted(list(set([s  for e in edlp+emlp # get the samples from the 1st
                                       for s in e.valsPerVar.itervalues().next().keys()])))
            print indent+'samples : ',samples
            colors = colors if colors else dict(zip(samples,
                                                    [r.kBlack,r.kRed,r.kBlue,r.kViolet,r.kGreen]))
            markers = markers if markers else dict(zip(samples, [24, 25, 26, 27, 28]))
            hprefix = pType.replace(' ','_')+'_'+l
            hDs = buildHistos(edlp, regions, hprefix+'_d')
            hMs = buildHistos(emlp, regions, hprefix+'_m')
            hRs = buildRatioHistos(hDs, hMs)
            for vt in hDs.keys() :
                can = plotHistos(histos1=hDs[vt], histos2=hMs[vt],
                                 histosRatio=hRs[vt], scaleXlabelsBot=1.5,
                                 label1='Davide', label2='Matt', canPrefix=hprefix,
                                 colors=colors, markers=markers, lineStyle1=1, lineStyle2=2)
                can.SaveAs(ioDir+hprefix+'_'+vt+'.png')
                print ioDir+hprefix+'_'+vt+'.png'
