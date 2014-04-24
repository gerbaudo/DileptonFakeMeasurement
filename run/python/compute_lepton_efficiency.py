#!/bin/env python

# compute eff(T|L) from MC truth
#
# davide.gerbaudo@gmail.com
# April 2014

import array
import collections
import glob
import math
import multiprocessing
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
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptTitle(1)
rootUtils.importRootCorePackages()
from datasets import datasets, setSameGroupForAllData
from SampleUtils import (colors
                         ,fastSamplesFromFilenames
                         ,guessSampleFromFilename
                         ,isBkgSample
                         ,markers)
from kin import computeMt

leptonSources = {'heavy' : 0, # see FakeLeptonSources.h
                 'light' : 1,
                 'conv' : 2,
                 'real' : 3,
                 'qcd' : 4,
                 'unknown' : 5
                 }
usage="""
Example usage:
%prog \\
 --verbose  \\
 --tag ${TAG} \\
 --output-dir ./out/fakerate/merged/el_eff_conv_${TAG}
 --lepton el
 --mode conv
 >& log/fakerate/el_eff_conv_${TAG}.log
"""
def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='out/fakerate/merged/lepton_efficiency')
    parser.add_option('-l', '--lepton', help='el or mu')
    parser.add_option('-m', '--mode', help='one of real,conv,hflf,hf,lf')
    parser.add_option('-t', '--tag', help='tag used to select the input files (e.g. Apr_04)')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    (options, args) = parser.parse_args()
    inputDir     = options.input_dir
    outputDir    = options.output_dir
    lepton       = options.lepton
    mode         = options.mode
    tag          = options.tag
    doFillHistos = options.fill_histos
    verbose      = options.verbose
    if not tag : parser.error('tag is a required option')
    validModes = ['real', 'conv', 'hflf', 'hf', 'lf']
    validLeptons = ['el', 'mu']
    if mode not in validModes : parser.error("invalid mode '%s', must be one of %s"%(mode, ','.join(validModes)))
    if lepton not in validLeptons : parser.error("invalid lepton %s, must be one of %s"%(lepton, ','.join(validLeptons)))

    tupleStems = {'real' : 'mcreal',
                  'conv' : 'mcconv',
                  'hflf' : 'mcqcd',
                  'hf'   : 'mcqcd',
                  'lf'   : 'mcqcd',
                  }
    treeNames = {'real' : 'RealExtractionRegion',
                 'conv' : 'ConversionExtractionRegion',
                 'hflf' : 'HfLfExtractionRegion',
                 'hf'   : 'HfLfExtractionRegion',
                 'lf'   : 'HfLfExtractionRegion'
                 }

    templateInputFilename = "*_%(m)s_tuple_%(t)s.root" % {'l':lepton, 'm':tupleStems[mode], 't':tag}
    templateOutputFilename =  "%(m)s_%(l)s_scale_histos.root" % {'l':lepton, 'm':mode}
    treeName = treeNames[mode]
    outputFileName = os.path.join(outputDir, templateOutputFilename)
    cacheFileName = outputFileName.replace('.root', '_cache.root')
    doFillHistos  = doFillHistos or not os.path.exists(cacheFileName)
    optionsToPrint = ['inputDir', 'outputDir', 'cacheFileName', 'lepton', 'mode', 'tag', 'doFillHistos']
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

    vars = ['pt', 'eta', 'pt_eta']
    groups = [g for g in samplesPerGroup.keys() if g!='data']
    #fill histos
    if doFillHistos :
        if verbose : print 'filling histograms'
        histosPerGroup = bookHistos(variables=vars, samples=groups, leptonType=lepton, mode=mode)
        for group in groups:
            filenames = filenamesPerGroup[group]
            histos = histosPerGroup[group]
            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            print "%s : %d entries"%(group, chain.GetEntries())
            fillHistos(chain, histos, lepton, mode, verbose, maxNevents=10000)
        writeHistos(cacheFileName, histosPerGroup, verbose)
    if verbose : print 'computing efficiencies'
    histosPerGroup = fetchHistos(cacheFileName, histoNames(vars, groups, lepton, mode), verbose)
    plotStackedHistos(histosPerGroup, outputDir, mode, verbose)
    efficiencies = computeEfficiencies(histosPerGroup)

    plotEfficiencies1d(efficiencies, outputDir, lepton, mode, verbose=verbose)
    plotEfficiencies2d(efficiencies, outputDir, lepton, mode, verbose=verbose)
    writeHistos(outputFileName, efficiencies, verbose)
    if verbose : print "saved efficiencies to %s" % outputFileName

#___________________________________________________

leptonTypes = ['tight', 'loose', 'real_tight', 'real_loose', 'fake_tight', 'fake_loose']

def histoname_electron_sf_vs_eta() : return 'sf_el_vs_eta'
def histoname_electron_sf_vs_pt() : return 'sf_el_vs_pt'

def fillHistos(chain, histos, lepton, mode, verbose=False, maxNevents=-1):
    nElecLoose, nElecTight = 0, 0
    totWeightLoose, totWeightTight = 0.0, 0.0
    validLepSources = {'real' : [leptonSources['real']],
                       'conv' : [leptonSources['conv']],
                       'hflf' : [leptonSources['light'], leptonSources['heavy']],
                       'hf'   : [leptonSources['heavy']],
                       'lf'   : [leptonSources['light']]
                       }[mode]
    def isRightSource(l) : return l.source in validLepSources
    def isRightLeptonType(l) : return l.isEl if lepton=='el' else l.isMu if lepton=='mu' else False
    for iEvent, event in enumerate(chain) :
        if maxNevents!=-1 and iEvent>maxNevents : break
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        l0, l1 = event.l0, event.l1
        def fill(l):
            if isRightLeptonType(l) and isRightSource(l) :
                l4m = r.TLorentzVector()
                l4m.SetPxPyPzE(l.px, l.py, l.pz, l.E)
                pt = l4m.Pt()
                eta = abs(l4m.Eta())
                histos['pt'    ]['loose'].Fill(pt, weight)
                histos['eta'   ]['loose'].Fill(eta, weight)
                histos['pt_eta']['loose'].Fill(pt, eta, weight)
                if l.isTight:
                    histos['pt'    ]['tight'].Fill(pt, weight)
                    histos['eta'   ]['tight'].Fill(eta, weight)
                    histos['pt_eta']['tight'].Fill(pt, eta, weight)
        fill(l0)
        fill(l1)

def histoName(var='', sample='', lType='', lSource='', mode='', looseOrTight=''):
    return 'h_'+sample+'_'+lType+'_'+lSource+'_'+mode+'_'+var+'_'+looseOrTight
def bookHistos(variables=[], samples=[], leptonType='', mode='') :
    "book a dict of histograms with keys [sample][var][loose,tight]"
    def histo(variable, hname):
        h = None
        ptBinEdges = np.array([10.0, 20.0, 35.0, 100.0])
        etaBinEdges = np.array([0.0, 1.37, 2.50])
        if   v=='pt'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif v=='eta'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        elif v=='pt_eta' : h = r.TH2F(hname, ';p_{T,l1}; #eta_{l1}',           len(ptBinEdges)-1,  ptBinEdges, len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    lt = leptonType
    return dict([(s,
                  dict([(v,
                         {'loose' : histo(variable=v, hname=histoName(v, s, lt, mode, 'loose')),
                          'tight' : histo(variable=v, hname=histoName(v, s, lt, mode, 'tight'))
                          })
                        for v in variables]))
                 for s in samples])
def histoNames(variables=[], samples=[], leptonType='', mode='') :
    def extractName(dictOrHist):
        "input must be either a dict or something with 'GetName'"
        isDict = type(dictOrHist) is dict
        return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()
    return extractName(bookHistos(variables, samples, leptonType, mode))

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
    return histos
def plotStackedHistos(histosPerGroup={}, outputDir='', mode='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    looseOrTight = first(first(histosPerGroup)).keys()
    histosPerName = dict([(mode+'_'+var+'_'+lt, # one canvas for each histo, so key with histoname w/out group
                           dict([(g, histosPerGroup[g][var][lt]) for g in groups]))
                          for var in variables for lt in looseOrTight])
    histosPerName = dict((k,v) for k, v in histosPerName.iteritems() # drop 2d, cannot draw stack
                         if not first(v).ClassName().startswith('TH2'))
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
        can.UseCurrentStyle()
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
        yMin, yMax = getMinMax([h for h in [totBkg, err_band] if h is not None])
        pm.SetMinimum(0.0)
        pm.SetMaximum(1.1*yMax)
        can.Update()
        topRightLabel(can, histoname, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]
        can.Update()
        can.SaveAs(os.path.join(outputDir, histoname+'.png'))

def computeEfficiencies(histosPerGroup={}) :
    groups = histosPerGroup.keys()
    vars = first(histosPerGroup).keys()
    num = dict((g, dict((v, histosPerGroup[g][v]['tight']) for v in vars)) for g in groups)
    den = dict((g, dict((v, histosPerGroup[g][v]['loose']) for v in vars)) for g in groups)
    eff = dict((g, dict((v, h.Clone(h.GetName().replace('tight', 'tight_over_loose')))
                         for v, h in num[g].iteritems()))
                for g in groups)
    [eff[g][v].Divide(den[g][v]) for g in groups for v in vars]
    # todo: set axes labels
    return eff

def plotEfficiencies1d(efficiencies, outputDir, lepton, mode, zoomIn=False, verbose=False):
    groups = efficiencies.keys()
    variables = [v for v in first(efficiencies).keys() if v not in ['pt_eta']]
    efficienciesPerVar = dict((v, dict((g, efficiencies[g][v]) for g in groups)) for v in variables)
    for var in variables:
        effs = efficienciesPerVar[var]
        padMaster = first(effs)
        framename = padMaster.GetName()
        def dropGroupName(s) :
            for g in groups : s = s.replace(g, '')
        dropGroupName(framename)
        can = r.TCanvas('c_'+framename, framename, 800, 600)
        can.cd()
        padMaster.SetStats(False)
        padMaster.Draw('axis')
        for g,h in effs.iteritems() :
            h.SetLineColor(colors[g] if g in colors else r.kBlack)
            h.SetMarkerColor(h.GetLineColor())
            h.SetMarkerStyle(markers[g] if g in markers else r.kFullCircle)
            h.Draw('ep same')
        minY, maxY = getMinMax(effs.values()) if zoomIn else (0.0, 1.0)
        padMaster.GetYaxis().SetRangeUser(min([0.0, minY]), 1.1*maxY)
        padMaster.SetMinimum(0.0)
        padMaster.SetMaximum(1.1*maxY)
        if maxY < 0.5 and zoomIn : padMaster.SetMaximum(0.25)
        padMaster.SetStats(False)
        xTitle = {'pt': 'p_{T} [GeV]', 'eta': '#eta'}[var]
        yTitle = '#varepsilon(T|L)'
        padMaster.SetTitle(';'+xTitle+';'+yTitle)
        frameTitle = "#varepsilon(T|L) %s %s"%(lepton, mode)
        drawLegendWithDictKeys(can, effs)
        topRightLabel(can, frameTitle, xpos=0.125, align=11)
        can.Update()
        can.SaveAs(os.path.join(outputDir, framename+'.png'))

def plotEfficiencies2d(efficiencies, outputDir, lepton, mode, verbose=False):
    groups = efficiencies.keys()
    variables = ['pt_eta']
    for g in groups:
        for v in variables:
            eff = efficiencies[g][v]
            framename = eff.GetName()
            can = r.TCanvas('c_'+framename, framename, 800, 600)
            can.cd()
            eff.SetStats(False)
            eff.SetTitle("#varepsilon(T|L) %s %s"%(lepton, mode))
            eff.Draw('colz')
            eff.Draw('text e same')
            can.Update()
            can.SaveAs(os.path.join(outputDir, framename+'.png'))

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
