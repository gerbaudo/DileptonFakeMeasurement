#!/bin/env python


import array
import collections
import glob
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
rootUtils.importRootCorePackages()
from datasets import datasets, setSameGroupForAllData
from SampleUtils import (colors
                         ,fastSamplesFromFilenames
                         ,guessSampleFromFilename
                         ,isBkgSample
                         ,isDataSample)
from kin import computeMt

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
    parser.add_option('-m', '--mode', default='hflf', help='either hflf or conv')
    parser.add_option('-t', '--tag', help='tag used to select the input files (e.g. Apr_04)')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    (options, args) = parser.parse_args()
    inputDir  = options.input_dir
    outputDir = options.output_dir
    lepton    = options.lepton
    mode      = options.mode
    tag       = options.tag
    verbose   = options.verbose
    if not tag : parser.error('tag is a required option')
    if lepton not in ['el', 'mu'] : parser.error("invalid lepton '%s'"%lepton)
    modes = ['hflf', 'conv'] if lepton=='el' else ['hflf']

    for mode in modes:
        print "running mode ",mode
        isConversion = mode=='conv'
        templateInputFilename = "*_%(mode)s_tuple_%(tag)s.root" % {'tag':tag, 'mode':mode}
        templateOutputFilename =  "%(mode)s_%(l)s_scale_histos.root" % {'mode':mode, 'l':lepton}
        treeName = 'HeavyFlavorControlRegion' if mode=='hflf' else 'ConversionControlRegion'
        outputFileName = os.path.join(outputDir, templateOutputFilename)
        cacheFileName = outputFileName.replace('.root', '_'+mode+'_cache.root')
        doFillHistograms = options.fill_histos or not os.path.exists(outputFileName)
        optionsToPrint = ['inputDir', 'outputDir', 'mode', 'tag', 'doFillHistograms']
        if verbose : print "options:\n"+'\n'.join(["%s : %s"%(o, eval(o)) for o in optionsToPrint])
        # collect inputs
        tupleFilenames = glob.glob(os.path.join(inputDir, templateInputFilename))
        samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
        samplesPerGroup = collections.defaultdict(list)
        filenamesPerGroup = collections.defaultdict(list)
        mkdirIfNeeded(outputDir)
        for s, f in zip(samples, tupleFilenames) :
            samplesPerGroup[s.group].append(s)
            filenamesPerGroup[s.group].append(f)
        vars = ['mt0', 'mt1', 'pt1', 'eta1']
        groups = samplesPerGroup.keys()
        #fill histos
        if doFillHistograms :
            histosPerGroup = bookHistos(vars, groups, mode=mode)
            histosPerSource = bookHistosPerSource(vars, leptonSources, mode=mode)
            for group in groups:
                isData = isDataSample(group)
                filenames = filenamesPerGroup[group]
                histosThisGroup = histosPerGroup[group]
                chain = r.TChain(treeName)
                [chain.Add(fn) for fn in filenames]
                print "%s : %d entries"%(group, chain.GetEntries())
                fillHistos(chain, histosThisGroup, histosPerSource, isConversion, isData, lepton, verbose)
            writeHistos(cacheFileName, histosPerGroup, histosPerSource, verbose)
        # compute scale factors
        histosPerGroup = fetchHistos(cacheFileName, histoNames(vars, groups, mode), verbose)
        histosPerSource = fetchHistos(cacheFileName, histoNamesPerSource(vars, leptonSources, mode), verbose)
        plotStackedHistos(histosPerGroup, outputDir, mode, verbose)
        plotStackedHistosSources(histosPerSource, outputDir, mode, verbose)
        sf_el_eta = subtractRealAndComputeScaleFactor(histosPerGroup, 'eta1', histoname_electron_sf_vs_eta(), verbose)
        sf_el_pt  = subtractRealAndComputeScaleFactor(histosPerGroup, 'pt1',  histoname_electron_sf_vs_pt(),  verbose)
        outputFile = r.TFile.Open(outputFileName, 'recreate')
        outputFile.cd()
        sf_el_eta.Write()
        sf_el_pt.Write()
        outputFile.Close()
        if verbose : print "saved scale factors to %s" % outputFileName

#___________________________________________________

leptonTypes = ['tight', 'loose', 'real_tight', 'real_loose', 'fake_tight', 'fake_loose']
allLeptonSources = ['heavy',   'light', 'conv',  'real',  'qcd', 'unknown'] # see FakeLeptonSources.h
leptonSources = [s for s in allLeptonSources if s not in ['qcd']] # qcd is just hf+lf
colorsSources = dict(zip(leptonSources, [r.kBlue-10, r.kMagenta-10, r.kRed-8, r.kGreen-6, r.kCyan-6, r.kGray+1]))

def enum2source(l): return allLeptonSources[l.source]
def histoname_electron_sf_vs_eta() : return 'sf_el_vs_eta'
def histoname_electron_sf_vs_pt() : return 'sf_el_vs_pt'

def fillHistos(chain, histosThisGroup, histosPerSource, isConversion, isData, lepton, verbose=False):
    nLoose, nTight = 0, 0
    totWeightLoose, totWeightTight = 0.0, 0.0
    for event in chain :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        tag, probe, met = event.l0, event.l1, event.met
        isSameSign = tag.charge*probe.charge > 0.
        isRightLep = probe.isEl if lepton=='el' else probe.isMu
        isTight = probe.isTight
        probeSource = probe.source
        sourceHf, sourceLf, sourceConv, sourceReal, sourceUnknown = 0, 1, 2, 3, 5 # see FakeLeptonSources.h
        isReal = probeSource==sourceReal and not isData
        isFake = not isReal and not isData
        def isBjet(j, mv1_80=0.3511) : return j.mv1 > mv1_80 # see SusyDefs.h
        # jets = event.jets # compute only if necessary
        # hasBjets = any(isBjet(j) for j in jets)
        tag4m, probe4m, met4m = r.TLorentzVector(), r.TLorentzVector(), r.TLorentzVector()
        tag4m.SetPxPyPzE(tag.px, tag.py, tag.pz, tag.E)
        probe4m.SetPxPyPzE(probe.px, probe.py, probe.pz, probe.E)
        met4m.SetPxPyPzE(met.px, met.py, met.pz, met.E)
        pt = probe4m.Pt()
        eta = abs(probe4m.Eta())
        mt0 = computeMt(tag4m, met4m)
        mt1 = computeMt(probe4m, met4m)
        isLowMt = mt1 < 40.0
        if (isSameSign or isConversion) and isRightLep and isLowMt:
            def incrementCounts(counts, weightedCounts):
                counts +=1
                weightedCounts += weight
            incrementCounts(nLoose, totWeightLoose)
            def fillHistosBySource(histos, probe):
                leptonSource = enum2source(probe)
                histos['mt1' ][leptonSource].Fill(mt1, weight)
                histos['pt1' ][leptonSource].Fill(pt,  weight)
                histos['eta1'][leptonSource].Fill(eta, weight)
            if not isData : fillHistosBySource(histosPerSource, probe)

            if isTight: incrementCounts(nTight, totWeightTight)
            histosThisGroup['mt0']['loose'].Fill(mt0, weight)
            histosThisGroup['mt1']['loose'].Fill(mt1, weight)
            def fill(lepType=''):
                histosThisGroup['pt1' ][lepType].Fill(pt, weight)
                histosThisGroup['eta1'][lepType].Fill(eta, weight)
            fill('loose')
            if isTight : fill('tight')
            if isReal : fill('real_loose')
            if isFake : fill('fake_loose')
            if isReal and isTight : fill('real_tight')
            if isFake and isTight : fill('fake_tight')
    if verbose:
        counterNames = ['nLoose', 'nTight', 'totWeightLoose', 'totWeightTight']
        print ', '.join(["%s : %.1f"%(c, eval(c)) for c in counterNames])

def histoNamePerSample(var, sample, leptonType, mode) : return 'h_'+var+'_'+sample+'_'+leptonType+'_'+mode
def histoNamePerSource(var, leptonSource, mode) : return 'h_'+var+'_'+leptonSource+'_'+mode
def bookHistos(variables, samples, leptonTypes=leptonTypes, mode='') :
    "book a dict of histograms with keys [sample][var][tight, loose, real_tight, real_loose]"
    def histo(variable, hname):
        h = None
        mtBinEdges = np.array([0.0, 20.0, 40.0, 60.0, 100.0, 200.0])
        ptBinEdges = np.array([10.0, 20.0, 35.0, 100.0])
        etaBinEdges = np.array([0.0, 1.37, 2.50])
        if   v=='mt0'     : h = r.TH1F(hname, ';m_{T}(tag,MET) [GeV]; entries/bin',   len(mtBinEdges)-1,  mtBinEdges)
        elif v=='mt1'     : h = r.TH1F(hname, ';m_{T}(probe,MET) [GeV]; entries/bin', len(mtBinEdges)-1,  mtBinEdges)
        elif v=='pt1'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif v=='eta1'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(s,
                  dict([(v,
                         dict([(lt, histo(variable=v, hname=histoNamePerSample(v, s, lt, mode)))
                               for lt in leptonTypes]))
                        for v in variables]))
                 for s in samples])
def bookHistosPerSource(variables, sources, mode=''):
    "book a dict of histograms with keys [var][lepton_source]"
    def histo(variable, hname):
        h = None
        mtBinEdges = np.array([0.0, 20.0, 40.0, 60.0, 100.0, 200.0])
        ptBinEdges = np.array([10.0, 20.0, 35.0, 100.0])
        etaBinEdges = np.array([0.0, 1.37, 2.50])
        if   variable=='mt0'     : h = r.TH1F(hname, ';m_{T}(tag,MET) [GeV]; entries/bin', len(mtBinEdges)-1,  mtBinEdges)
        elif variable=='mt1'     : h = r.TH1F(hname, ';m_{T}(probe,MET) [GeV]; entries/bin', len(mtBinEdges)-1,  mtBinEdges)
        elif variable=='pt1'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif variable=='eta1'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(v,
                  dict([(s, histo(variable=v, hname=histoNamePerSource(v, s, mode)))
                        for s in leptonSources]))
                 for v in variables])
def extractName(dictOrHist):
    "input must be either a dict or something with 'GetName'"
    isDict = type(dictOrHist) is dict
    return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()
def histoNames(variables, samples, mode) :
    return extractName(bookHistos(variables, samples, mode=mode))
def histoNamesPerSource(variables, samples, mode) :
    return extractName(bookHistosPerSource(variables, leptonSources, mode=mode))


def writeHistos(outputFileName='', histosPerGroup={}, histosPerSource={}, verbose=False):
    outputFile = r.TFile.Open(outputFileName, 'recreate')
    outputFile.cd()
    if verbose : print "writing to %s"%outputFile.GetName()
    def write(dictOrObj):
        isDict = type(dictOrObj) is dict
        if isDict:
            for v in dictOrObj.values():
                write(v)
        else:
            dictOrObj.Write()
    write(histosPerGroup)
    write(histosPerSource)
    outputFile.Close()
def fetchHistos(fileName='', histoNames={}, verbose=False):
    "given a dict of input histonames, return the same dict, but with histo instead of histoname"
    inputFile = r.TFile.Open(fileName)
    if verbose : print "fetching histograms from %s"%inputFile.GetName()
    def fetch(dictOrName):
        isDict = type(dictOrName) is dict
        return dict([(k, fetch(v)) for k,v in dictOrName.iteritems()]) if isDict else inputFile.Get(dictOrName)
    histos = fetch(histoNames)
    return histos
def plotStackedHistos(histosPerGroup={}, outputDir='', mode='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    leptonTypes = first(first(histosPerGroup)).keys()
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
        yMin, yMax = getMinMax([h for h in [totBkg, data, err_band] if h is not None])
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
    for var, histos in histosPerVar.iteritems():
        canvasBasename = mode+'_region_'+var
        missingSources = [s for s, h in histosPerVar.iteritems() if not h]
        if missingSources:
            if verbose : print "skip %s, missing histos for %s"%(histoname, str(missingSources))
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
            print 'adding ',s
            h.SetFillColor(colorsSources[s] if s in colorsSources else r.kOrange)
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
def subtractRealAndComputeScaleFactor(histosPerGroup={}, variable='', outhistoname='', verbose=False):
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
    return ratio

def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
    print "computeStatErr2 use the one in rootUtils"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}

def buildErrBandGraph(histo_tot_bkg, err2s) :
    print "buildErrBandGraph use the one in rootUtils"
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
