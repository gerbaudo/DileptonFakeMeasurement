#!/bin/env python


import array
import collections
import glob
import math
import numpy as np
import optparse
import os
import pprint
from utils import (first
                   ,mkdirIfNeeded
                   )
import rootUtils
from rootUtils import (drawLegendWithDictKeys
                       ,getBinContents
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
                         ,isBkgSample)
from kin import computeMt

def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='./out/fake_el_scale_factor')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False)
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    (options, args) = parser.parse_args()
    inputDir = options.input_dir
    outputDir = options.output_dir
    tag = 'Apr_02'
    outputFileName = os.path.join(outputDir, "fake_el_scale_histos_%s.root"%tag)
    doFillHistograms = options.fill_histos or not os.path.exists(outputFileName)
    verbose = options.verbose
    mkdirIfNeeded(outputDir)
    tupleFilenames = glob.glob(os.path.join(inputDir, '*_fake_tuple_'+tag+'.root'))
    samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
    samplesPerGroup = collections.defaultdict(list)
    filenamesPerGroup = collections.defaultdict(list)
    for s, f in zip(samples, tupleFilenames) :
        samplesPerGroup[s.group].append(s)
        filenamesPerGroup[s.group].append(f)
    vars = ['pt1', 'eta1']
    groups = samplesPerGroup.keys()
    if doFillHistograms :
        histosPerGroup = bookHistos(vars, groups)
        for group in groups:
            treeName = 'HeavyFlavorControlRegion'
            filenames = filenamesPerGroup[group]
            histos = histosPerGroup[group]
            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            print "%s : %d entries"%(group, chain.GetEntries())
            fillHistos(chain, histos, verbose)
        writeHistos(outputFileName, histosPerGroup, verbose)
    histosPerGroup = fetchHistos(outputFileName, histoNames(vars, groups), verbose)
    plotStackedHistos(histosPerGroup, outputDir, verbose)
    subtractRealAndComputeScaleFactor(histosPerGroup, 'eta1', verbose)
    subtractRealAndComputeScaleFactor(histosPerGroup, 'pt1', verbose)


def fillHistos(chain, histos, verbose=False):
    totNelec, totWeightLoose, totWeightTight = 0, 0.0, 0.0
    nElecTight = 0
    nOutRange, wOutRange = 0, 0.0
    for event in chain :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        tag, probe, met = event.l0, event.l1, event.met
        isSameSign = tag.charge*probe.charge > 0.
        isEl, isTight = probe.isEl, probe.isTight
        isReal = probe.source==3 # see FakeLeptonSources.h
        probe4m, met4m = r.TLorentzVector(), r.TLorentzVector()
        probe4m.SetPxPyPzE(probe.px, probe.py, probe.pz, probe.E)
        met4m.SetPxPyPzE(met.px, met.py, met.pz, met.E)
        pt = probe4m.Pt()
        eta = abs(probe4m.Eta())
        mt = computeMt(probe4m, met4m)
        if not isSameSign : continue
        if not isEl : continue
        if mt > 100.0 : print "mt ",mt
        if mt > 40.0 : continue
        if pt<10.0 or pt>100.0 :
            nOutRange += 1
            wOutRange += weight
        totNelec += 1
        totWeightLoose += weight
        def fill(lepType=''):
            histos['pt1'][lepType].Fill(pt, weight)
            histos['eta1'][lepType].Fill(eta, weight)
        fill('loose')
        if isTight:
            fill('tight')
            totWeightTight += weight
            nElecTight += 1
        if isReal : fill('real_loose')
        if isTight and isReal : fill('real_tight')

    if verbose:
        counterNames = ['totNelec', 'nElecTight', 'totWeightLoose', 'totWeightTight', 'nOutRange', 'wOutRange']
        print ', '.join(["%s : %.1f"%(c, eval(c)) for c in counterNames])

def histoName(var, sample, leptonType) : return 'h_'+var+'_'+sample+'_'+leptonType
def bookHistos(variables, samples) :
    "book a dict of histograms with keys [sample][var][tight, loose, real_tight, real_loose]"
    def histo(variable, hname):
        h = None
        ptBinEdges = np.array([10.0, 20.0, 35.0, 100.0])
        etaBinEdges = np.array([0.0, 1.37, 2.50])
        if   v=='pt1'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif v=='eta1'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    lepTypes = ['tight', 'loose', 'real_tight', 'real_loose']
    return dict([(s,
                  dict([(v,
                         dict([(lt, histo(variable=v, hname=histoName(v, s, lt)))
                               for lt in lepTypes]))
                        for v in variables]))
                 for s in samples])
def histoNames(variables, samples) :
    def extractName(dictOrHist):
        "input must be either a dict or something with 'GetName'"
        isDict = type(dictOrHist) is dict
        return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()
    return extractName(bookHistos(variables, samples))
def writeHistos(outputFileName='', histosPerGroup={}, verbose=False):
    outputFile = r.TFile.Open(outputFileName, 'recreate')
    outputFile.cd()
    if verbose : print "writing to %s"%outputFile.GetName()
    def write(dictOrObj):
        isDict = type(dictOrObj) is dict
        if isDict:
            for v in dictOrObj.values():
                write(v)
        else:
            if verbose : print dictOrObj.GetName()
            dictOrObj.Write()
    write(histosPerGroup)
    outputFile.Close()
def fetchHistos(fileName='', histoNames={}, verbose=False):
    "given a dict of input histonames, return the same dict, but with histo instead of histoname"
    inputFile = r.TFile.Open(fileName)
    if verbose : print "fetching histograms from %s"%inputFile.GetName()
    def fetch(dictOrName):
        isDict = type(dictOrName) is dict
        return dict([(k, fetch(v)) for k,v in dictOrName.iteritems()]) if isDict else inputFile.Get(dictOrName)
    histos = fetch(histoNames)
    if verbose : print "fetched histos:\n%s"%pprint.pformat(histos)
    return histos
def plotStackedHistos(histosPerGroup={}, outputDir='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    leptonTypes = first(first(histosPerGroup)).keys()
    histosPerName = dict([(var+'_'+lt, # one canvas for each histo, so key with histoname w/out group
                           dict([(g, histosPerGroup[g][var][lt]) for g in groups]))
                          for var in variables for lt in leptonTypes])
    for histoname, histosPerGroup in histosPerName.iteritems():
        missingGroups = [g for g, h in histosPerGroup.iteritems() if not h]
        if missingGroups:
            if verbose : print "skip %s, missing histos for %s"%(histoname, str(missingGroups))
            continue
        bkgHistos = dict([(g, h) for g, h in histosPerGroup.iteritems() if isBkgSample(g)])
        totBkg = summedHisto(bkgHistos.values())
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
        data = histosPerGroup['data']
        if data and data.GetEntries():
            data.SetMarkerStyle(r.kFullDotLarge)
            data.Draw('p same')
        yMin, yMax = getMinMax([h for h in [totBkg, data] if h is not None])
        pm.SetMinimum(0.0)
        pm.SetMaximum(1.1*yMax)
        can.Update()
        topRightLabel(can, histoname, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, bkgHistos, opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]+[data]
        can.Update()
        can.SaveAs(os.path.join(outputDir, histoname+'.png'))
def subtractRealAndComputeScaleFactor(histosPerGroup={}, variable='', verbose=False):
    "efficiency scale factor"
    groups = histosPerGroup.keys()
    leptonTypes = ['tight', 'loose', 'real_tight', 'real_loose']
    histosPerType = dict([(lt,
                           dict([(g,
                                  histosPerGroup[g][variable][lt])
                                 for g in groups]))
                          for lt in leptonTypes])
    for lt in leptonTypes :
        histosPerType[lt]['totSimBkg'] = summedHisto([histo for group,histo in histosPerType[lt].iteritems() if isBkgSample(group)])

    dataTight = histosPerType['tight']['data'     ]
    simuTight = histosPerType['tight']['totSimBkg']
    dataLoose = histosPerType['loose']['data'     ]
    simuLoose = histosPerType['loose']['totSimBkg'] # todo: for mc, build fake_tight, fake_loose rather than subtracting (smaller error)
    dataTight.Add(histosPerType['real_tight']['totSimBkg'], -1.0)
    simuTight.Add(histosPerType['real_tight']['totSimBkg'], -1.0)
    dataLoose.Add(histosPerType['real_loose']['totSimBkg'], -1.0)
    simuLoose.Add(histosPerType['real_loose']['totSimBkg'], -1.0)
    dataTight.Divide(dataLoose)
    simuTight.Divide(simuLoose)
    print "eff(T|L) vs. ",variable
    print "efficiency data : ",getBinContents(dataTight)
    print "efficiency simu : ",getBinContents(simuTight)
    dataTight.Divide(simuTight)
    print "scale factor data/simu: ",getBinContents(dataTight)




if __name__=='__main__':
    main()
