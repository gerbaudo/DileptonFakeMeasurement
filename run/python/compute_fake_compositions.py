#!/bin/env python

# script to compute the fake lepton compositions for the weighted average

# This script computes the fake lepton composition from the
# 'fake-ntuples' produced by MeasureFakeRate2.
#
# input: out/fakerate/*_ssinc1j_tuple_*.root
# output: tbd
#
# davide.gerbaudo@gmail.com
# April 2014

import array
import collections
import glob
import kin
import math
import numpy as np
import optparse
import os
import pprint
from utils import (dictSum
                   ,first
                   ,mkdirIfNeeded
                   )
import rootUtils
from rootUtils import (drawLegendWithDictKeys
                       ,getBinContents
                       ,getBinErrors
                       ,getMinMax
                       ,importRoot
                       ,importRootCorePackages
                       ,summedHisto
                       ,topRightLabel)
r = rootUtils.importRoot()
r.gROOT.SetStyle('Plain')
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
rootUtils.importRootCorePackages()
from datasets import datasets, setSameGroupForAllData
from SampleUtils import (fastSamplesFromFilenames
                         ,guessSampleFromFilename
                         ,isBkgSample
                         ,isDataSample)
import SampleUtils
from kin import computeMt
import fakeUtils as fakeu

usage="""
Example usage:
%prog \\
 --verbose  \\
 --tag ${TAG} \\
 --output-dir ./out/fakerate/el_sf_${TAG}
 >& log/fakerate/el_sf_${TAG}.log
"""
def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='./out/fake_el_scale_factor', help='dir for plots')
    parser.add_option('-l', '--lepton', default='el', help='either el or mu')
    parser.add_option('-t', '--tag', help='tag used to select the input files (e.g. Apr_04)')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    (options, args) = parser.parse_args()
    inputDir  = options.input_dir
    outputDir = options.output_dir
    lepton    = options.lepton
    tag       = options.tag
    verbose   = options.verbose
    if not tag : parser.error('tag is a required option')
    if lepton not in ['el', 'mu'] : parser.error("invalid lepton '%s'"%lepton)

    templateInputFilename = "*_ssinc1j_tuple_%(tag)s.root" % {'tag':tag}
    templateOutputFilename =  "%(l)s_composition_histos.root" % {'l':lepton}
    treeName = 'SameSign1jetControlRegion'
    outputFileName = os.path.join(outputDir, templateOutputFilename)
    cacheFileName = outputFileName.replace('.root', '_cache.root')
    doFillHistograms = options.fill_histos or not os.path.exists(outputFileName)
    optionsToPrint = ['inputDir', 'outputDir', 'tag', 'doFillHistograms']
    if verbose : print "options:\n"+'\n'.join(["%s : %s"%(o, eval(o)) for o in optionsToPrint])
    # collect inputs
    print '----> input files ',os.path.join(inputDir, templateInputFilename)
    tupleFilenames = glob.glob(os.path.join(inputDir, templateInputFilename))
    samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
    samplesPerGroup = collections.defaultdict(list)
    filenamesPerGroup = collections.defaultdict(list)
    mkdirIfNeeded(outputDir)
    for s, f in zip(samples, tupleFilenames) :
        samplesPerGroup[s.group].append(s)
        filenamesPerGroup[s.group].append(f)
    vars = ['pt', 'eta']
    groups = samplesPerGroup.keys()
    #fill histos
    if doFillHistograms :
        histosPerGroupPerSource = bookHistosPerSamplePerSource(vars, groups, leptonSources)
        for group in groups:
            isData = isDataSample(group)
            filenames = filenamesPerGroup[group]
            histosThisGroupPerSource = histosPerGroupPerSource[group]
            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            print "%s : %d entries"%(group, chain.GetEntries())
            fillHistos(chain, histosThisGroupPerSource, isData, lepton, group, verbose)
        writeHistos(cacheFileName, histosPerGroupPerSource, verbose)
#     # compute scale factors
#     histosPerGroup = fetchHistos(cacheFileName, histoNames(vars, groups, mode), verbose)
#     histosPerSource = fetchHistos(cacheFileName, histoNamesPerSource(vars, leptonSources, mode), verbose)
#     histosPerSamplePerSource = fetchHistos(cacheFileName, histoNamesPerSamplePerSource(vars, groups, leptonSources, mode), verbose)
#     plotStackedHistos(histosPerGroup, outputDir, mode, verbose)
#     plotStackedHistosSources(histosPerSource, outputDir, mode, verbose)
#     plotPerSourceEff(histosPerVar=histosPerSource, outputDir=outputDir, lepton=lepton, mode=mode, verbose=verbose)
#     for g in groups:
#         hps = dict((v, histosPerSamplePerSource[v][g])for v in vars)
#         plotPerSourceEff(histosPerVar=hps, outputDir=outputDir, lepton=lepton, mode=mode, sample=g, verbose=verbose)
#     sf_el_eta = subtractRealAndComputeScaleFactor(histosPerGroup, 'eta1', histoname_electron_sf_vs_eta(), outputDir, mode, verbose)
#     sf_el_pt  = subtractRealAndComputeScaleFactor(histosPerGroup, 'pt1',  histoname_electron_sf_vs_pt(),  outputDir, mode, verbose)
#     outputFile = r.TFile.Open(outputFileName, 'recreate')
#     outputFile.cd()
#     sf_el_eta.Write()
#     sf_el_pt.Write()
#     outputFile.Close()
#     if verbose : print "saved scale factors to %s" % outputFileName

#___________________________________________________

leptonTypes = fakeu.leptonTypes()
allLeptonSources = fakeu.allLeptonSources()
leptonSources = fakeu.leptonSources()
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def allSelections() :
    return ['ssinc1j']+[srcr+'_'+ll+'_'+nj
                        for srcr in ['sr','cr'] for ll in ['ee','em','mm'] for nj in ['eq1j','ge2j']]
def histoname_electron_sf_vs_eta() : return 'sf_el_vs_eta'
def histoname_electron_sf_vs_pt() : return 'sf_el_vs_pt'

def getSelection(l0, l1, jets, met):
    def isBjet(j, mv1_80=0.3511) : return j.mv1 > mv1_80 # see SusyDefs.h
    def isForwardJet(j) : return j.p4.Pt()>30.0 and abs(j.p4.Eta())<2.4 # should be detEta but not stored
    fabs = math.fabs
    clJets = [j for j in jets if not isBjet(j) and not isForwardJet(j) for j in jets]
    hasFbJets = len(clJets)!=len(jets)
    nClJets = len(clJets)
    sel = None
    if nClJets==0:
        return sel
    l0pt    = l0.p4.Pt()
    l1pt    = l1.p4.Pt()
    j0      = clJets[0]
    mll     = (l0.p4 + l1.p4).M()
    ht      = kin.computeHt(met.p4, [l0.p4, l1.p4]+[j.p4 for j in clJets])
    metrel  = kin.computeMetRel(met.p4, [l0.p4, l1.p4]+[j.p4 for j in clJets])
    mtl0    = kin.computeMt(l0.p4, met.p4)
    mtl1    = kin.computeMt(l1.p4, met.p4)
    mtmin   = min([mtl0, mtl1])
    mtmax   = max([mtl0, mtl1])
    mlj     = kin.computeMlj(l0.p4, l1.p4, j0.p4)
    detall  = fabs(l0.p4.Eta() - l1.p4.Eta())
    mljj = kin.computeMljj(l0.p4, l1.p4, clJets[0].p4, clJets[1].p4) if nClJets>1 else None
    ll = kin.getDilepType(l0, l1)
    nj = 'eq1j' if nClJets==1 else 'ge2j'
    if (ll=='mm' and nj=='eq1j'
        and l0pt    >  30.0
        and l1pt    >  20.0
        and detall  <   1.5
        and mtmax   > 100.0
        and ht      > 200.0):
        sel = 'sr_mm_eq1j' if mlj<90.0 else 'cr_mm_eq1j'
    elif (ll=='mm' and nj=='ge2j'
          and l0pt    >  30.0
          and l1pt    >  30.0
          and detall  <   1.5
          and ht      > 220.0):
        sel ='sr_mm_ge2j' if mljj<120.0 else 'cr_mm_ge2j'
    elif (ll=='em' and nj=='eq1j'
          and l0pt    >  30.0
          and l1pt    >  30.0
          and detall  <   1.5
          and ht      > 200.0
          and mtmax   > 110.0):
        sel = 'sr_ee_eq1j' if mlj<90.0 else 'cr_ee_eq1j'
    elif (ll=='em' and nj=='ge2j'
          and l0pt    >  30.0
          and l1pt    >  30.0
          and detall  <   1.5
          and ht      > 200.0
          and mtmax   > 110.0):
        sel = 'sr_em_ge2j' if mljj<120.0 else 'cr_em_ge2j'
    elif (ll=='ee' and nj=='eq1j'
          and l0pt    >  30.0
          and l1pt    >  20.0
          and fabs(mll-91.2) > 10.0
          and metrel  >  55.0
          and ht      > 200.0):
        sel = 'sr_ee_eq1j' if mlj<90.0 else 'cr_ee_eq1j'
    elif (ll=='ee' and nj=='ge2j'
          and l0pt    >  30.0
          and l1pt    >  20.0
          and fabs(mll-91.2) > 10.0
          and metrel  >  30.0
          and mtmax   > 100.0):
        sel = 'sr_ee_ge2j' if mljj<120.0 else 'cr_ee_ge2j'
    return sel



def fillHistos(chain, histosThisGroupPerSource, isData, lepton, group, verbose=False):
    "expect histos[group][sel][source][var][loose,tight]"
    normFactor = 3.2 if group=='heavyflavor' else 1.0 # bb/cc hand-waving normalization factor, see notes 2014-04-17
    nLepFilled = 0
    for event in chain :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        weight = weight*normFactor
        l0, l1 = kin.addTlv(event.l0), kin.addTlv(event.l1)
        met = kin.addTlv(event.met)
        jets = [kin.addTlv(j) for j in event.jets]
        sourceReal = 3 # see FakeLeptonSources.h
        l0IsFake = l0.source!=sourceReal and not isData
        l1IsFake = l0.source!=sourceReal and not isData
        atLeastOneIsFake = l0IsFake or l1IsFake
        if not atLeastOneIsFake : continue
        selection = getSelection(l0, l1, jets, met)
        def fillHistosBySource(lep):
            isTight = lep.isTight
            source = lep.source
            leptonSource = enum2source(lep)
            isReal = source==sourceReal and not isData
            isFake = not isReal and not isData
            sourceIsKnown = not isData
            isRightLep = lep.isMu and lepton=='mu' or lep.isEl and lepton=='el'
            def fill():
                pt, eta = lep.p4.Pt(), abs(lep.p4.Eta())
                histosThisGroupPerSource['ssinc1j'][leptonSource]['pt' ]['loose'].Fill(pt,  weight)
                histosThisGroupPerSource['ssinc1j'][leptonSource]['eta']['loose'].Fill(eta, weight)
                if selection:
                    histosThisGroupPerSource[selection][leptonSource]['pt' ]['loose'].Fill(pt,  weight)
                    histosThisGroupPerSource[selection][leptonSource]['eta']['loose'].Fill(eta, weight)
            filled = False
            if isRightLep and sourceIsKnown and isFake :
                fill()
                filled = True
            return filled
        if fillHistosBySource(l0) : nLepFilled += 1
        if fillHistosBySource(l1) : nLepFilled += 1
    if verbose:
        print "filled histos for %d leptons"%nLepFilled

def histoNamePerSamplePerSource(var, group, sel, source, tightOrLoose):
    return 'h_'+var+'_'+group+'_'+sel+'_'+source+'_'+tightOrLoose
def bookHistosPerSamplePerSource(variables, samples, sources):
    "book a dict of histograms with keys [group][sel][source][var][tight, loose]"
    def histo(variable, hname):
        h = None
        ptBinEdges = fakeu.ptBinEdges()
        etaBinEdges = fakeu.etaBinEdges()
        if   variable=='pt'  : h = r.TH1F(hname, ';p_{T,l} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif variable=='eta' : h = r.TH1F(hname, ';#eta_{l}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    print 'variables ',variables
    print 'samples ',samples
    print 'sources ',sources
    return dict([(g,
                  dict([(sel,
                         dict([(s,
                                dict([(v,
                                       {'loose' : histo(variable=v, hname=histoNamePerSamplePerSource(v, g, sel, s, 'loose')),
                                        'tight' : histo(variable=v, hname=histoNamePerSamplePerSource(v, g, sel, s, 'tight')),
                                        })
                                      for v in variables]))
                               for s in leptonSources]))
                        for sel in allSelections()]))
                 for g in samples])

def extractName(dictOrHist):
    "input must be either a dict or something with 'GetName'"
    isDict = type(dictOrHist) is dict
    return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()
def histoNames(variables, samples, mode) :
    return extractName(bookHistos(variables, samples, mode=mode))
def histoNamesPerSource(variables, samples, mode) :
    return extractName(bookHistosPerSource(variables, leptonSources, mode=mode))
def histoNamesPerSamplePerSource(variables, samples, leptonSources, mode) :
    return extractName(bookHistosPerSamplePerSource(variables, samples, leptonSources, mode=mode))

def writeHistos(outputFileName='', histosPerSamplePerSource={}, verbose=False):
    rootUtils.writeObjectsToFile(outputFileName, histosPerSamplePerSource, verbose)
def fetchHistos(fileName='', histoNames={}, verbose=False):
    return rootUtils.fetchObjectsFromFile(fileName, histoNames, verbose)


def plotStackedHistos(histosPerGroup={}, outputDir='', mode='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    leptonTypes = first(first(histosPerGroup)).keys()
    colors = SampleUtils.colors
    histosPerName = dict([(mode+'_'+var+'_'+lt, # one canvas for each histo, so key with histoname w/out group
                           dict([(g, histosPerGroup[g][var][lt]) for g in groups]))
                          for var in variables for lt in leptonTypes])
    for histoname, histosPerGroup in histosPerName.iteritems():
        missingGroups = [g for g, h in histosPerGroup.iteritems() if not h]
        if missingGroups:
            if verbose : print "skip %s, missing histos for %s"%(histoname, str(missingGroups))
            continue
        bkgHistos = dict([(g, h) for g, h in histosPerGroup.iteritems() if isBkgSample(g)])
        totBkg = summedHisto(bkgHistos.values())
        err_band = buildErrBandGraph(totBkg, computeStatErr2(totBkg))
        emptyBkg = totBkg.Integral()==0
        if emptyBkg:
            if verbose : print "empty backgrounds, skip %s"%histoname
            continue
        can = r.TCanvas('c_'+histoname, histoname, 800, 600)
        can.cd()
        pm = totBkg # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update() # necessary to fool root's dumb object ownership
        stack = r.THStack('stack_'+histoname,'')
        can.Update()
        r.SetOwnership(stack, False)
        for s, h in bkgHistos.iteritems() :
            h.SetFillColor(colors[s] if s in colors else r.kOrange)
            h.SetDrawOption('bar')
            h.SetDirectory(0)
            stack.Add(h)
        stack.Draw('hist same')
        err_band.Draw('E2 same')
        data = histosPerGroup['data']
        if data and data.GetEntries():
            data.SetMarkerStyle(r.kFullDotLarge)
            data.Draw('p same')
        yMin, yMax = getMinMax([h for h in [totBkg, data, err_band] if h])
        pm.SetMinimum(0.0)
        pm.SetMaximum(1.1*yMax)
        can.Update()
        topRightLabel(can, histoname, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]+[data]
        can.Update()
        can.SaveAs(os.path.join(outputDir, histoname+'.png'))
def plotStackedHistosSources(histosPerVar={}, outputDir='', mode='', verbose=False):
    variables = histosPerVar.keys()
    sources = first(histosPerVar).keys()
    colors = colorsFillSources
    for var in variables:
        histos = dict((s, histosPerVar[var][s]['loose']) for s in sources)
        canvasBasename = mode+'_region_'+var+'_loose'
        missingSources = [s for s, h in histos.iteritems() if not h]
        if missingSources:
            if verbose : print "skip %s, missing histos for %s"%(var, str(missingSources))
            continue
        totBkg = summedHisto(histos.values())
        err_band = buildErrBandGraph(totBkg, computeStatErr2(totBkg))
        emptyBkg = totBkg.Integral()==0
        if emptyBkg:
            if verbose : print "empty backgrounds, skip %s"%canvasBasename
            continue
        can = r.TCanvas('c_'+canvasBasename, canvasBasename, 800, 600)
        can.cd()
        pm = totBkg # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update() # necessary to fool root's dumb object ownership
        stack = r.THStack('stack_'+canvasBasename,'')
        can.Update()
        r.SetOwnership(stack, False)
        for s, h in histos.iteritems() :
            h.SetFillColor(colors[s] if s in colors else r.kOrange)
            h.SetDrawOption('bar')
            h.SetDirectory(0)
            stack.Add(h)
        stack.Draw('hist same')
        err_band.Draw('E2 same')
        yMin, yMax = getMinMax([h for h in [totBkg, err_band] if h is not None])
        pm.SetMinimum(0.0)
        pm.SetMaximum(1.1*yMax)
        can.Update()
        topRightLabel(can, canvasBasename, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(histos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]
        can.Update()
        can.SaveAs(os.path.join(outputDir, canvasBasename+'.png'))

def plotPerSourceEff(histosPerVar={}, outputDir='', lepton='', mode='', sample='', verbose=False, zoomIn=True):
    "plot efficiency for each source (and 'anysource') as a function of each var; expect histos[var][source][loose,tight]"
    variables = histosPerVar.keys()
    sources = [s for s in first(histosPerVar).keys() if s!='real'] # only fake eff really need a scale factor
    colors = colorsLineSources
    for var in filter(lambda x : x in ['pt1', 'eta1'], histosPerVar.keys()):
        histosPerSource = dict((s, histosPerVar[var][s]) for s in sources)
        canvasBasename = mode+'_efficiency_'+lepton+'_'+var+("_%s"%sample if sample else '')
        missingSources = [s for s, h in histosPerSource.iteritems() if not h['loose'] or not h['tight']]
        if missingSources:
            if verbose : print "skip %s, missing histos for %s"%(var, str(missingSources))
            continue
        anySourceLoose = summedHisto([h['loose'] for h in histosPerSource.values()])
        anySourceTight = summedHisto([h['tight'] for h in histosPerSource.values()])
        anySourceLoose.SetName(histoNamePerSource(var, 'any', 'loose', mode))
        anySourceTight.SetName(histoNamePerSource(var, 'any', 'tight', mode))
        histosPerSource['any'] = { 'loose' : anySourceLoose, 'tight' : anySourceTight }
        emptyBkg = anySourceLoose.Integral()==0 or anySourceTight.Integral()==0
        if emptyBkg:
            if verbose : print "empty backgrounds, skip %s"%canvasBasename
            continue
        def computeEfficiencies(histosPerSource={}) :
            sources = histosPerSource.keys()
            num = dict((s, histosPerSource[s]['tight']) for s in sources)
            den = dict((s, histosPerSource[s]['loose']) for s in sources)
            eff = dict((s, h.Clone(h.GetName().replace('tight', 'tight_over_loose')))
                       for s, h in num.iteritems())
            [eff[s].Divide(den[s]) for s in sources]
            return eff

        effs = computeEfficiencies(histosPerSource)
        can = r.TCanvas('c_'+canvasBasename, canvasBasename, 800, 600)
        can.cd()
        pm = first(effs) # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update()
        for s, h in effs.iteritems() :
            h.SetMarkerColor(colors[s] if s in colors else r.kBlack)
            h.SetLineColor(h.GetMarkerColor())
            h.SetLineWidth(2*h.GetLineWidth())
            h.SetMarkerStyle(markersSources[s] if s in markersSources else r.kDot)
            h.Draw('ep same')
            h.SetDirectory(0)
        #pprint.pprint(effs)
        yMin, yMax = getMinMax(effs.values())
        pm.SetMinimum(0.0)
        pm.SetMaximum(0.25 if yMax < 0.5 and zoomIn else 1.1)
        can.Update()
        topRightLabel(can, canvasBasename, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, effs, opt='lp')
        can.RedrawAxis()
        can._histos = effs
        can.Update()
        can.SaveAs(os.path.join(outputDir, canvasBasename+'.png'))

def subtractRealAndComputeScaleFactor(histosPerGroup={}, variable='', outhistoname='', outputDir='./', mode='', verbose=False):
    "efficiency scale factor"
    groups = histosPerGroup.keys()
    histosPerType = dict([(lt,
                           dict([(g,
                                  histosPerGroup[g][variable][lt])
                                 for g in groups]))
                          for lt in leptonTypes])
    for lt in leptonTypes :
        histosPerType[lt]['totSimBkg'] = summedHisto([histo for group,histo in histosPerType[lt].iteritems() if isBkgSample(group)])

    simuTight = histosPerType['fake_tight']['totSimBkg']
    simuLoose = histosPerType['fake_loose']['totSimBkg']
    dataTight = histosPerType['tight'     ]['data'     ]
    dataLoose = histosPerType['loose'     ]['data'     ]
    # subtract real contribution from data
    # _Note to self_: currently estimating the real contr from MC; in
    # the past also used iterative corr, which might be more
    # appropriate in cases like here, where the normalization is
    # so-so.  Todo: investigate the normalization.
    dataTight.Add(histosPerType['real_tight']['totSimBkg'], -1.0)
    dataLoose.Add(histosPerType['real_loose']['totSimBkg'], -1.0)
    dataTight.Divide(dataLoose)
    simuTight.Divide(simuLoose)
    print "eff(T|L) vs. ",variable
    def formatFloat(floats): return ["%.4f"%f for f in floats]
    print "efficiency data : ",formatFloat(getBinContents(dataTight))
    print "efficiency simu : ",formatFloat(getBinContents(simuTight))
    ratio = dataTight.Clone(outhistoname)
    ratio.SetDirectory(0)
    ratio.Divide(simuTight)
    print "sf    data/simu : ",formatFloat(getBinContents(ratio))
    print "            +/- : ",formatFloat(getBinErrors(ratio))
    can = r.TCanvas('c_'+outhistoname, outhistoname, 800, 600)
    botPad, topPad = rootUtils.buildBotTopPads(can)
    can.cd()
    topPad.Draw()
    topPad.cd()
    pm = dataTight
    pm.SetStats(0)
    pm.Draw('axis')
    xAx, yAx = pm.GetXaxis(), pm.GetYaxis()
    xAx.SetTitle('')
    xAx.SetLabelSize(0)
    yAx.SetRangeUser(0.0, 0.25)
    textScaleUp = 1.0/topPad.GetHNDC()
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetTitle('#epsilon(T|L)')
    yAx.SetTitleOffset(yAx.GetTitleOffset()/textScaleUp)
    simuTight.SetLineColor(r.kRed)
    simuTight.SetMarkerStyle(r.kOpenCross)
    simuTight.SetMarkerColor(simuTight.GetLineColor())
    dataTight.Draw('same')
    simuTight.Draw('same')
    leg = drawLegendWithDictKeys(topPad, {'data':dataTight, 'simulation':simuTight}, legWidth=0.4)
    leg.SetHeader('scale factor '+mode+' '+('electron' if '_el_'in outhistoname else 'muon' if '_mu_' in outhistoname else ''))
    can.cd()
    botPad.Draw()
    botPad.cd()
    ratio.SetStats(0)
    ratio.Draw()
    textScaleUp = 1.0/botPad.GetHNDC()
    xAx, yAx = ratio.GetXaxis(), ratio.GetYaxis()
    yAx.SetRangeUser(0.0, 2.0)
    xAx.SetTitle({'pt1':'p_{T}', 'eta1':'|#eta|'}[variable])
    yAx.SetNdivisions(-202)
    yAx.SetTitle('Data/Sim')
    yAx.CenterTitle()
    xAx.SetLabelSize(textScaleUp*0.04)
    xAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    refLine = rootUtils.referenceLine(xAx.GetXmin(), xAx.GetXmax())
    refLine.Draw()
    can.Update()
    can.SaveAs(os.path.join(outputDir, mode+'_'+outhistoname+'.png'))
    can.SaveAs(os.path.join(outputDir, mode+'_'+outhistoname+'.eps'))
    return ratio

def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
#     print "computeStatErr2 use the one in rootUtils"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}

def buildErrBandGraph(histo_tot_bkg, err2s) :
#     print "buildErrBandGraph use the one in rootUtils"
    h = histo_tot_bkg
    bins = range(1, 1+h.GetNbinsX())
    x = np.array([h.GetBinCenter (b) for b in bins])
    y = np.array([h.GetBinContent(b) for b in bins])
    ex_lo = ex_hi = np.array([0.5*h.GetBinWidth(b) for b in bins])
    ey_lo, ey_hi = np.sqrt(err2s['down']), np.sqrt(err2s['up'])
    gr = r.TGraphAsymmErrors(len(bins), x, y, ex_lo, ex_hi, ey_lo, ey_hi)
    gr.SetMarkerSize(0)
    gr.SetFillStyle(3004)
    gr.SetFillColor(r.kGray+3)
    gr.SetLineWidth(2)
    return gr


if __name__=='__main__':
    main()
