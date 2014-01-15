#!/bin/env python

# Optimize the final selection for [ee,em,mm] x [1j, 2j]
#
# Input: tuples from susy::wh::TupleMaker
#
# davide.gerbaudo@gmail.com
# Jan 2014

import datetime
import glob
import math
fabs = math.fabs
import optparse
import os
from rootUtils import importRoot, buildBotTopPads, summedHisto, binContentsWithUoflow, cloneAndFillHisto, cumEffHisto, maxSepVerticalLine
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

from kin import (phi_mpi_pi,
                 addTlv,
                 computeMt, computeHt, computeMetRel,
                 getDilepType,
                 computeMt2, computeMt2j, computeMljj,
                 thirdLepZcandidateIsInWindow)
from utils import (getCommandOutput,
                   guessLatestTagFromLatestRootFiles,
                   guessMonthDayTagFromLastRootFile,
                   isMonthDayTag,
                   dictSum,
                   first,
                   rmIfExists,
                   linearTransform,
                   cumsum,
                   mergeOuter
                   )
from SampleUtils import isSigSample, colors

def optimizeSelection() :
    inputdir, options = parseOptions()
    tag = pickTag(inputdir, options)
    sigFiles, bkgFiles = getInputFilenames(inputdir, tag, options)
    allSamples = dictSum(sigFiles, bkgFiles)
    vars = variablesToPlot()
    histos = bookHistos(vars, allSamples.keys(), options.ll, options.nj)
    fillHistos(histos, dictSum(sigFiles, bkgFiles), options.ll, options.nj)
    bkgHistos = dict((s, h) for s, h in histos.iteritems() if s in bkgFiles.keys())
    sigHistos = dict((s, h) for s, h in histos.iteritems() if s in sigFiles.keys())
    plotHistos(bkgHistos, sigHistos)
    printSummary(histos)
    
def parseOptions() :
    usage="""%prog [options] dir
    Example:
    %prog -v out/susysel
    """
    lls, njs = ['ee','em','mm'], ['eq1j','ge2j']
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('--ll', help='dilepton (default all)')
    parser.add_option('--nj', help='njet multiplicity (default all)')
    parser.add_option("-s", "--sample-regexp", dest="samples", default='.*', help="consider only matching samples (default '.*')")
    parser.add_option("-e", "--exclude-regexp", dest="exclude", default=None, help="exclude matching samples")
    parser.add_option('-t', '--tag', help='production tag; by default the latest one')
    parser.add_option('-v', '--verbose', action='store_true', help='print details')
    parser.add_option('-d', "--debug", action='store_true', help='print even more details')
    (options, args) = parser.parse_args()
    if len(args) != 1 : parser.error("incorrect number of arguments")
    def validateMultiOpt(o, opts, defaults) :
        v = getattr(opts, o) if hasattr(opts, o) else None
        if v and v not in defaults : parser.error("%s must be in %s (got '%s')"%(o, str(defaults), v))
        setattr(opts, o, defaults if v is None else [v])
    validateMultiOpt('ll', options, lls)
    validateMultiOpt('nj', options, njs)
    inputdir = args[0]
    return inputdir, options

def pickTag(inputdir, options) :
    tag = options.tag if options.tag else guessLatestTagFromLatestRootFiles(inputdir, options.debug)
    tag = tag if tag else guessMonthDayTagFromLastRootFile(inputdir, options.debug) # can get it wrong if from single file
    tag = tag.strip('_') # leading/trailing separators are not part of the tag
    if not isMonthDayTag(tag) : print "warning, non-standard tag might lead to bugs '%s'"%tag
    if options.verbose : print "using tag %s"%tag
    return tag

def getInputFilenames(inputdir, tag, options) :
    sig, bkg = dict(), dict()
    for f in glob.glob(inputdir+'/*'+tag+'.root') :
        sample = os.path.splitext(os.path.basename(f))[0]
        coll = sig if isSigSample(sample) else bkg
        assert sample not in coll,"%s already found\n%s\n%s"%(sample, f, coll[sample])
        coll[sample] = f
    if options.verbose : print 'input files:\n',sig,'\n',bkg
    assert sig and bkg, "missing signals or backgrounds"
    return sig, bkg

def variablesToPlot() :
    return ['pt0','pt1'] #,'mll','mtmin',]#'mtmax','mtllmet','ht','metrel','dphill','detall',
            #'mt2j','mljj','dphijj','detajj']

def llnjKey(ll, nj) : return "%s_%s"%(ll, nj)
def histoSuffix(sample, ll, nj) : return "%s_%s_%s"%(sample, ll, nj)

def bookHistos(variables, samples, lls, njs) :
    "book a dict of histograms with keys [sample][ll_nj][var]"
    def histo(variable, suffix) :
        s = suffix
        twopi = +2.0*math.pi
        if   v=='pt0'     : return r.TH1F('h_pt0_'    +s, ';p_{T,l0} [GeV]; entries/bin',          25, 0.0, 250.0)
        elif v=='pt1'     : return r.TH1F('h_pt1_'    +s, ';p_{T,l1} [GeV]; entries/bin',          25, 0.0, 250.0)
        elif v=='mll'     : return r.TH1F('h_mll_'    +s, ';m_{l0,l1} [GeV]; entries/bin',         25, 0.0, 250.0)
        elif v=='mtmin'   : return r.TH1F('h_mtmin_'  +s, ';m_{T,min}(l, MET) [GeV]; entries/bin', 25, 0.0, 250.0)
        elif v=='mtmax'   : return r.TH1F('h_mtmax_'  +s, ';m_{T,max}(l, MET) [GeV]; entries/bin', 25, 0.0, 250.0)
        elif v=='mtllmet' : return r.TH1F('h_mtllmet_'+s, ';m_{T}(l+l, MET) [GeV]; entries/bin',   25, 0.0, 250.0)
        elif v=='ht'      : return r.TH1F('h_ht_'     +s, ';H_{T} [GeV]; entries/bin',             25, 0.0, 250.0)
        elif v=='metrel'  : return r.TH1F('h_metrel_' +s, ';MET_{rel} [GeV]; entries/bin',         25, 0.0, 250.0)
        elif v=='dphill'  : return r.TH1F('h_dphill_' +s, ';#Delta#phi(l, l) [rad]; entries/bin',  25, 0.0, twopi)
        elif v=='detall'  : return r.TH1F('h_detall_' +s, ';#Delta#eta(l, l); entries/bin',        25, 0.0, +3.0 )
        elif v=='mt2j'    : return r.TH1F('h_mt2j_'   +s, ';m^{J}_{T2} [GeV]; entries/bin',        25, 0.0, 500.0)
        elif v=='mljj'    : return r.TH1F('h_mljj_'   +s, ';m_{ljj} [GeV]; entries/bin',           25, 0.0, 500.0)
        elif v=='dphijj'  : return r.TH1F('h_dphijj_' +s, ';#Delta#phi(j, j); entries/bin',        25, 0.0, twopi)
        elif v=='detajj'  : return r.TH1F('h_detajj_' +s, '#Delta#eta(j, j); entries/bin',         25, 0.0, +3.0 )
        else : print "unknown variable %s"%v
    return dict([(s,
                  dict([(llnjKey(ll, nj),
                         dict([(v, histo(v, histoSuffix(s, ll, nj))) for v in variables]))
                         for ll in lls for nj in njs]))
                 for s in samples])

def fillHistos(histos, files, lls, njs) :
    treename = 'SusySel'
    for sample, filename in files.iteritems() :
        histosSample = histos[sample]
        tree = r.TFile.Open(filename).Get(treename)
        print "processing %s (%d entries) %s"%(sample, tree.GetEntries(), datetime.datetime.now())
        for event in tree :
            l0, l1, met, pars = addTlv(event.l0), addTlv(event.l1), addTlv(event.met), event.pars
            jets, lepts = [addTlv(j) for j in event.jets], [addTlv(l) for l in event.lepts]
            ll = getDilepType(l0, l1)
            nJets = len(jets)
            nj = 'eq1j' if nJets==1 else 'ge2j'
            if ll not in lls or nj not in njs : continue
            pt0 = l0.p4.Pt()
            pt1 = l1.p4.Pt()
            mll  = (l0.p4 + l1.p4).M()
            mtllmet = computeMt(l0.p4 + l1.p4, met.p4)
            ht      = computeHt(met.p4, [l0.p4, l1.p4]+[j.p4 for j in jets])
            metrel  = computeMetRel(met.p4, [l0.p4, l1.p4]+[j.p4 for j in jets])
            mtl0    = computeMt(l0.p4, met.p4)
            mtl1    = computeMt(l1.p4, met.p4)
            mtmin   = min([mtl0, mtl1])
            mtmax   = max([mtl0, mtl1])
            dphill  = abs(phi_mpi_pi(l0.p4.DeltaPhi(l1.p4))) 
            detall  = fabs(l0.p4.Eta() - l1.p4.Eta())
            l3Veto  =  not thirdLepZcandidateIsInWindow(l0, l1, lepts)
            if nJets >1 :
                j0, j1 = jets[0], jets[1]
                mt2j   = computeMt2j(l0.p4, l1.p4, j0.p4, j1.p4, met.p4)
                mljj   = computeMljj(l0.p4, l1.p4, j0.p4, j1.p4)
                dphijj = fabs(phi_mpi_pi(j0.p4.DeltaPhi(j1.p4)))
                detajj = fabs(j0.p4.Eta() - j1.p4.Eta())
            if passSelection(pt0, pt1, mll, mtllmet, met, jets, lepts, ll, nj) :
                varHistos = histosSample[llnjKey(ll, nj)]
                varValues = dict([(v, eval(v)) for v in variablesToPlot()])
                fillVarHistos(varHistos, varValues, pars.weight, nj)
                
def passSelection(l0pt, l1pt, mll, mtllmet, ht, metrel, l3Veto, ll, nj) :
    def passLepPT(l0pt, l1pt, ll) :
        return l0pt > 30.0 and l1pt > (20.0 if ll=='mm' else 0.0)
    def passZveto(mll, ll) :
        return True if ll!='ee' else fabs(mll - 91.2) > 10.0
    def passMtLlMetMin(mtllmet) :
        return mtllmet > (150.0 if ll=='ee' else 140.0 if ll=='em' else 100 if ll=='mm' else 0.0)
    def passHtMin(ht) : return  ht > 200.0
    def passMetRelMin(metrel, ll) : return metrel > (50.0 if ll!='mm' else 0.0)
    return (passLepPT(l0pt, l1pt, ll)
            and passZveto(mll, ll)
            and passMtLlMetMin(mtllmet)
            and passHtMin(ht)
            and passMetRelMin(metrel, ll)
            and l3Veto
            )

def fillVarHistos(varHistos, varValues, weight, nj) :
    exclVars = ['mt2j','mljj','dphijj','detajj'] if nj < 2 else []
    vars = [v for v in varHistos.keys() if v not in exclVars]
    for v in vars :
            varHistos[v].Fill(varValues[v], weight)

def plotHistos(bkgHistos, sigHistos) :
    llnjs = first      (sigHistos).keys()
    vars  = first(first(sigHistos)).keys()
    for llnj in llnjs :
        for var in vars :
            plotVar(dict((s, bkgHistos[s][llnj][var]) for s in bkgHistos.keys()),
                    dict((s, sigHistos[s][llnj][var]) for s in sigHistos.keys()),
                    llnj+'_'+var)

def plotVar(bkgHistos, sigHistos, llnjvar) :
    def preferredSignal(signals):
        pref = 'Herwigpp_sM_wA_noslep_notauhad_WH_2Lep_1'
        return pref if pref in signals else first(sorted(signals))
    signalSample = preferredSignal(sigHistos.keys())
    allHistos = bkgHistos.values() + [sigHistos[signalSample],]
    allHistosEmpty = all([h.GetEntries()==0 for h in allHistos])
    if allHistosEmpty : return
    can = r.TCanvas('can_'+llnjvar, llnjvar, 800, 800)
    botPad, topPad = buildBotTopPads(can, splitFraction=0.75)
    totBkg = summedHisto(bkgHistos.values())
    can._totBkg = totBkg
    can.cd()
    botPad.Draw()
    drawBottom(botPad, totBkg, bkgHistos, sigHistos[signalSample], llnjvar)
    can.cd()
    topPad.Draw()
    drawTop(topPad, totBkg, sigHistos[signalSample])
    can.Update()
    outFilename = llnjvar+'.png'
    rmIfExists(outFilename) # avoid root warnings
    can.SaveAs(outFilename)

def drawBottom(pad, totBkg, bkgHistos, sigHisto, llnjvar) :
    pad.cd()
    totBkg.Draw('axis')
    stack = r.THStack('stack_'+llnjvar,'')
    for s, h in bkgHistos.iteritems() :
        h.SetFillColor(colors[s] if s in colors else r.kOrange)
        h.SetDrawOption('bar')
        stack.Add(h)
    stack.Draw('hist same')
    hTot = stack.GetHistogram()
    sigHisto.SetLineColor(r.kRed)
    sigHisto.SetLineWidth(2*sigHisto.GetLineWidth())
    # legend todo
    sigHisto.Draw('same')
    def writeLabel(can, label, font='') :
        tex = r.TLatex(0.0, 0.0, '')
        tex.SetNDC()
        if font : tex.SetTextFont(font)
        tex.SetTextAlign(31)
        tex.DrawLatex(1.0-pad.GetTopMargin(), 1.0-pad.GetRightMargin(), label)
        return tex
    pad._lab = writeLabel(llnjvar)
    pad._histos = [stack]
    pad.Update()

def drawTop(pad, hBkg, hSig) :
    nxS, nxB = hSig.GetNbinsX()+1, hBkg.GetNbinsX()+1  # TH1 starts from 1
    assert nxS==nxB,"histos with differen binning (%d!=%d)"%(nxS,nxB)
    pad.SetGridy()
    bcS, bcB = binContentsWithUoflow(hSig), binContentsWithUoflow(hBkg)
    leftToRight = True
    bcLS, bcLB = cumsum(mergeOuter(bcS), leftToRight), cumsum(mergeOuter(bcB), leftToRight),
    leftToRight = False    
    bcRS, bcRB = cumsum(mergeOuter(bcS), leftToRight), cumsum(mergeOuter(bcB), leftToRight)
    zn = r.RooStats.NumberCountingUtils.BinomialExpZ
    bkgUnc = 0.3
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
    pad._obj = [hCeS, hCeB, mark, gr, midline]

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

        

    
def printSummary(histos) :
    refVar = 'pt0'
    for sample, hs in histos.iteritems() :
        for sel, hss in hs.iteritems() :
            print "%s %s %f"%(sample, sel, hss[refVar].Integral())
            
    

if __name__=='__main__' :
    optimizeSelection()
