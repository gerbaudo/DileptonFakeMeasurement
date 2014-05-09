#!/bin/env python

# script to plot (from the fake ntuples) the variables (iso, pv) used to define the tight leptons
#
# davide.gerbaudo@gmail.com
# May 2014

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
                   ,rmIfExists
                   )
import rootUtils
from rootUtils import (drawLegendWithDictKeys
                       ,getBinContents
                       ,getBinErrors
                       ,getMinMax
                       ,importRoot
                       ,importRootCorePackages
                       ,summedHisto
                       ,topRightLabel
                       ,rightLegend)
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
import fakeUtils as fakeu
from compute_fake_el_scale_factor import histoname_electron_sf_vs_eta
from buildWeightedMatrix import fetchSfHistos
from compute_fake_el_scale_factor import computeStatErr2, buildErrBandGraph

usage="""
Example usage:
%prog \\
 --verbose  \\
 --tag ${TAG} \\
 --output-dir ./out/fakerate/el_sf_${TAG}
 >& log/fakerate/el_sf_${TAG}.log

 TODO
"""

def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='./out/tight_variables_plots', help='dir for plots')
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

    mode = 'zmmejets'
    # ./out/fakerate/periodJ.physics_Egamma.PhysCont_zmmejets_tuple_May_07.root
    templateInputFilename = "*_%(mode)s_tuple_%(tag)s.root" % {'tag':tag, 'mode':mode}
    templateOutputFilename =  "%(mode)s_%(l)s_tight_plots.root" % {'mode':mode, 'l':lepton}
    treeName = 'ZmmeVetoPlusJetsRegion' if mode=='zmmejets' else 'unknownTreeName'
    outputFileName = os.path.join(outputDir, templateOutputFilename)
    cacheFileName = outputFileName.replace('.root', '_'+mode+'_cache.root')
    doFillHistograms = options.fill_histos or not os.path.exists(cacheFileName)
    optionsToPrint = ['inputDir', 'outputDir', 'mode', 'tag', 'doFillHistograms']
    if verbose : print "options:\n"+'\n'.join(["%s : %s"%(o, eval(o)) for o in optionsToPrint])
    # collect inputs
    if verbose : print 'input files ',os.path.join(inputDir, templateInputFilename)
    tupleFilenames = glob.glob(os.path.join(inputDir, templateInputFilename))
    samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
    samplesPerGroup = collections.defaultdict(list)
    filenamesPerGroup = collections.defaultdict(list)
    mkdirIfNeeded(outputDir)
    for s, f in zip(samples, tupleFilenames) :
        samplesPerGroup[s.group].append(s)
        filenamesPerGroup[s.group].append(f)
    vars = ['pt','eta','d0sig','z0SinTheta','etCone','ptCone','etConeCorr','ptConeCorr']
    groups = samplesPerGroup.keys()
    #fill histos
    if doFillHistograms :
        histosPerGroup = bookHistosPerGroup(vars, groups)
        for group in groups:
            isData = isDataSample(group)
            filenames = filenamesPerGroup[group]
            histosThisGroup = histosPerGroup[group]
            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            print "%s : %d entries"%(group, chain.GetEntries())
            fillHistos(chain, histosThisGroup, isData, lepton, group, verbose)
        writeHistos(cacheFileName, histosPerGroup, verbose)
    # compute scale factors
    histosPerGroup = fetchHistos(cacheFileName, histoNames(vars, groups), verbose)
    plotStackedHistos(histosPerGroup, outputDir, mode, verbose)
#___________________________________________________

leptonTypes = fakeu.leptonTypes()
allLeptonSources = fakeu.allLeptonSources()
leptonSources = fakeu.leptonSources()
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def allSelections() :
    return ['incl']
     #[srcr+'_'+ll+'_'+nj for srcr in ['sr','cr'] for ll in ['ee','em','mm'] for nj in ['eq1j','ge2j']]

def writeHistos(outputFileName='', histosNames={}, verbose=False):
    rootUtils.writeObjectsToFile(outputFileName, histosNames, verbose)

def fetchHistos(fileName='', histoNames={}, verbose=False):
    return rootUtils.fetchObjectsFromFile(fileName, histoNames, verbose)

def histoNamePerSample(sample, var, tightOrLoose) : return "h_%s_%s_%s"%(sample, var, tightOrLoose)

def bookHistosPerGroup(variables, groups) :
    "book a dict of histograms with keys [sample][var][tight,loose]"
    def histo(variable, hname):
        h = None
        ptBinEdges = fakeu.ptBinEdges()
        etaBinEdges = fakeu.etaBinEdges()
        # you are here!!!!
        if   v=='pt'          : h = r.TH1F(hname, ';p_{T}(probe e) [GeV]; entries/bin',      len(ptBinEdges)-1, ptBinEdges)
        elif v=='eta'         : h = r.TH1F(hname, ';#eta(probe e); entries/bin',             len(etaBinEdges)-1, etaBinEdges)
        elif v=='d0sig'       : h = r.TH1F(hname, ';d_{0 sig}(probe e); entries/bin',        50, -50.0, +50.0)
        elif v=='z0SinTheta'  : h = r.TH1F(hname, ';z_{0}sin#theta(probe e); entries/bin',   50, -50.0, +50.0)
        elif v=='etCone'      : h = r.TH1F(hname, ';E_{T,cone}(probe e) [GeV]; entries/bin', 60, -10.0, +50.0)
        elif v=='ptCone'      : h = r.TH1F(hname, ';p_{T,cone}(probe e) [GeV]; entries/bin', 50,   0.0, +50.0)
        elif v=='etConeCorr'  : h = r.TH1F(hname, ';E_{T, cone, corr}(probe e); entries/bin',60, -10.0, +50.0)
        elif v=='ptConeCorr'  : h = r.TH1F(hname, ';p_{T, cone, corr}(probe e); entries/bin',50,   0.0, +50.0)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(g,
                  dict([(v,
                         {'loose' : histo(variable=v, hname=histoNamePerSample(g, v, 'loose')),
                          'tight' : histo(variable=v, hname=histoNamePerSample(g, v, 'tight'))
                          })
                        for v in variables]))
                 for g in groups])

def extractName(dictOrHist):
    "input must be either a dict or something with 'GetName'"
    isDict = type(dictOrHist) is dict
    return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()

def histoNames(variables, samples) :
    return extractName(bookHistosPerGroup(variables, samples))

def fillHistos(chain, histosThisGroup, isData, lepton, group, verbose=False):
    nLoose, nTight = 0, 0
    totWeightLoose, totWeightTight = 0.0, 0.0
    for event in chain :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        weight = weight
        probe = kin.addTlv(event.l1)
        met = kin.addTlv(event.met)
        isRightLep = probe.isEl if lepton=='el' else probe.isMu
        isTight = probe.isTight
        pt = probe.p4.Pt()
        eta = abs(probe.p4.Eta())
        d0Sig, z0Sin = probe.d0Signif, probe.z0SinTheta
        etCone, ptCone = probe.etCone, probe.ptCone
        etConeCorr, ptConeCorr = probe.etConeCorr, probe.ptConeCorr 
        selection = True
        if selection:
            def incrementCounts(counts, weightedCounts):
                counts +=1
                weightedCounts += weight
            def fill(tightOrLoose):
                histosThisGroup['pt'        ][tightOrLoose].Fill(pt, weight)
                histosThisGroup['eta'       ][tightOrLoose].Fill(eta, weight)
                histosThisGroup['d0sig'     ][tightOrLoose].Fill(d0Sig, weight)
                histosThisGroup['z0SinTheta'][tightOrLoose].Fill(z0Sin, weight)
                histosThisGroup['etCone'    ][tightOrLoose].Fill(etCone, weight)
                histosThisGroup['ptCone'    ][tightOrLoose].Fill(ptCone, weight)
                histosThisGroup['etConeCorr'][tightOrLoose].Fill(etConeCorr, weight)
                histosThisGroup['ptConeCorr'][tightOrLoose].Fill(ptConeCorr, weight)
            incrementCounts(nLoose, totWeightLoose)
            fill('loose')
            if isTight:
                incrementCounts(nTight, totWeightTight)
                fill('tight')
    if verbose:
        counterNames = ['nLoose', 'nTight', 'totWeightLoose', 'totWeightTight']
        print ', '.join(["%s : %.1f"%(c, eval(c)) for c in counterNames])

def plotStackedHistos(histosPerGroup={}, outputDir='', mode='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    colors = SampleUtils.colors
    histosPerName = dict([(var+'_'+tol, # one canvas for each histo, so key with histoname w/out group
                           dict([(g, histosPerGroup[g][var][tol]) for g in groups]))
                          for var in variables for tol in ['tight','loose']])
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
        topRightLabel(can, histoname+', '+mode, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]+[data]
        can.Update()
        can.SaveAs(os.path.join(outputDir, histoname+'.png'))


if __name__=='__main__':
    main()
