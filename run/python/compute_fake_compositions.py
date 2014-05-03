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
from kin import computeMt
import fakeUtils as fakeu
import plotParametrizedFractions

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
    doFillHistograms = options.fill_histos or not os.path.exists(cacheFileName)
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
    # compute and plot fractions
    histosPerGroupPerSource = fetchHistos(cacheFileName, histoNamesPerSamplePerSource(vars, groups, leptonSources))
    for sel in allSelections():
        for var in vars:
            hs, groups = histosPerGroupPerSource, histosPerGroupPerSource.keys()
            groups = [g for g in groups if g!='data' and g!='higgs']
            histosHeavy = dict((g, hs[g][sel]['heavy'][var]['loose']) for g in groups)
            histosLight = dict((g, hs[g][sel]['light'][var]['loose']) for g in groups)
            histosConv  = dict((g, hs[g][sel]['conv' ][var]['loose']) for g in groups)
            normalizeHistos   = plotParametrizedFractions.normalizeHistos
            plotStackedHistos = plotParametrizedFractions.plotStackedHistos

            frameTitle = 'hf elec: '+sel+' loose;'+var
            canvasName = 'elec_hf'+sel+'_'+var+'_den'
            plotStackedHistos(histosHeavy, canvasName, outputDir, frameTitle)

            frameTitle = 'lf elec: '+sel+' loose;'+var
            canvasName = 'elec_lf'+sel+'_'+var+'_den'
            plotStackedHistos(histosHeavy, canvasName, outputDir, frameTitle)

            frameTitle = 'conv elec: '+sel+' loose;'+var
            canvasName = 'elec_conv'+sel+'_'+var+'_den'
            plotStackedHistos(histosConv, canvasName, outputDir, frameTitle)

            # normalize and draw fractions (den only)
            histos = dict([(k+'_heavy',  h) for k,h in histosHeavy.iteritems()] +
                          [(k+'_light',  h) for k,h in histosLight.iteritems()] +
                          [(k+'_conv', h) for k,h in histosConv.iteritems()])
            normalizeHistos(histos)
            histos = {'heavy':histosHeavy, 'light':histosLight, 'conv':histosConv}
            frameTitle = 'elec: '+sel+';'+var
            canvasBaseName = 'elec_fake'+sel+'_'+var+'_frac'
            plotFractionsStacked(histos, canvasBaseName+'_stack', outputDir, frameTitle)


def plotFractionsStacked(histos={}, canvasName='', outputDir='./', frameTitle='title;p_{T} [GeV]; fraction'):
    "plot stack of histos[source][group]"
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    can.SetRightMargin(2.0*can.GetRightMargin())
    stack = r.THStack('stack_'+canvasName, '')
    leg = rightLegend(can)
    leg.SetBorderSize(0)
    colors = SampleUtils.colors
    sources = sorted(histos.keys())
    groups  = sorted(first(histos).keys())
    fillStyles = {'light':3365, 'heavy':3356, 'conv':3663}
    assert len(sources)<4,"can only plot up to three sources, %s"%str(sources)
    def setHistoAtts(histo, group, colors=colors, fillStyle=None) :
        h, g = histo, group
        h.SetMarkerSize(0)
        h.SetFillColor(colors[g])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
        if fillStyle : h.SetFillStyle(fillStyle)
    for s in sources:
        for g in groups:
            h = histos[s][g]
            setHistoAtts(h, g, fillStyle=fillStyles[s])
            stack.Add(h)
    for s in sources[::-1]: # stack goes b-t, legend goes t-b
        leg.AddEntry(r.TObject(), s, '')
        def niceLabel(group) : return 'bb/cc' if group=='heavyflavor' else group
        for g in groups[::-1] : leg.AddEntry(histos[s][g], niceLabel(g) , 'F')
    stack.Draw('hist')
    leg.Draw('same')
    tex = r.TLatex()
    tex.SetNDC(True)
    tex.DrawLatex(0.1, 0.9, 'fractions: '+frameTitle[:frameTitle.find(';')-1])
    can.Update() # force stack to create canMaster
    canMaster = stack.GetHistogram()
    canMaster.SetTitle(frameTitle)
    canMaster.Draw('axis same')
    can._graphical_objects = [stack, canMaster, leg, tex] + [h for h in stack.GetStack()]
    can.Update()
    for ext in ['png'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)

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


def getSelectionChannelOnly(l0, l1, jets, met):
    nClJets = len(jets) #note to self: these are already central-light jets b/c of passCommonCriteria
    sel = None
    ll = kin.getDilepType(l0, l1)
    nj = 'eq1j' if nClJets==1 else 'ge2j'
    mlj  = kin.computeMlj(l0.p4, l1.p4, jets[0].p4)
    mljj = kin.computeMljj(l0.p4, l1.p4, jets[0].p4, jets[1].p4) if nClJets>1 else None
    if (ll=='mm' and nj=='eq1j'):
        sel = 'sr_mm_eq1j' if mlj<90.0 else 'cr_mm_eq1j'
    elif (ll=='mm' and nj=='ge2j'):
        sel ='sr_mm_ge2j' if mljj<120.0 else 'cr_mm_ge2j'
    elif (ll=='em' and nj=='eq1j'):
        sel = 'sr_ee_eq1j' if mlj<90.0 else 'cr_ee_eq1j'
    elif (ll=='em' and nj=='ge2j'):
        sel = 'sr_em_ge2j' if mljj<120.0 else 'cr_em_ge2j'
    elif (ll=='ee' and nj=='eq1j'):
        sel = 'sr_ee_eq1j' if mlj<90.0 else 'cr_ee_eq1j'
    elif (ll=='ee' and nj=='ge2j'):
        sel = 'sr_ee_ge2j' if mljj<120.0 else 'cr_ee_ge2j'
    return sel

def getSelection(l0, l1, jets, met):
    nClJets = len(jets)
    sel = None
    l0pt    = l0.p4.Pt()
    l1pt    = l1.p4.Pt()
    j0      = jets[0]
    mll     = (l0.p4 + l1.p4).M()
    ht      = kin.computeHt(met.p4, [l0.p4, l1.p4]+[j.p4 for j in jets])
    metrel  = kin.computeMetRel(met.p4, [l0.p4, l1.p4]+[j.p4 for j in jets])
    mtl0    = kin.computeMt(l0.p4, met.p4)
    mtl1    = kin.computeMt(l1.p4, met.p4)
    mtmin   = min([mtl0, mtl1])
    mtmax   = max([mtl0, mtl1])
    mlj     = kin.computeMlj(l0.p4, l1.p4, j0.p4)
    detall  = fabs(l0.p4.Eta() - l1.p4.Eta())
    mljj = kin.computeMljj(l0.p4, l1.p4, jets[0].p4, jets[1].p4) if nClJets>1 else None
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
        selection = getSelectionChannelOnly(l0, l1, jets, met)
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
def histoNamesPerSamplePerSource(variables, samples, leptonSources) :
    return extractName(bookHistosPerSamplePerSource(variables, samples, leptonSources))

def writeHistos(outputFileName='', histosPerSamplePerSource={}, verbose=False):
    rootUtils.writeObjectsToFile(outputFileName, histosPerSamplePerSource, verbose)
def fetchHistos(fileName='', histoNames={}, verbose=False):
    return rootUtils.fetchObjectsFromFile(fileName, histoNames, verbose)


if __name__=='__main__':
    main()
