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
import utils
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
    parser.add_option('-r', '--region', help='one of the regions for which we saved the fake ntuples')
    parser.add_option('-t', '--tag', help='tag used to select the input files (e.g. Apr_04)')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    (options, args) = parser.parse_args()
    inputDir  = options.input_dir
    outputDir = options.output_dir
    lepton    = options.lepton
    region    = options.region
    tag       = options.tag
    verbose   = options.verbose
    if not tag : parser.error('tag is a required option')
    if lepton not in ['el', 'mu'] : parser.error("invalid lepton '%s'"%lepton)
    filestems, treenames = utils.verticalSlice(fakeu.tupleStemsAndNames)
    regions = filestems
    assert region in regions,"invalid region '%s', must be one of %s"%(region, str(regions))
    # you are here : implement different regions
    templateInputFilename = "*_%(region)s_tuple_%(tag)s.root" % {'tag':tag, 'region':region}
    templateOutputFilename =  "%(region)s_%(l)s_tight_plots.root" % {'region':region, 'l':lepton}
    treeName = treenames[regions.index(region)]
    outputDir = outputDir+'/'+region+'/'+lepton # split the output in subdirectories, so we don't overwrite things
    mkdirIfNeeded(outputDir)
    outputFileName = os.path.join(outputDir, templateOutputFilename)
    cacheFileName = outputFileName.replace('.root', '_'+region+'_cache.root')
    doFillHistograms = options.fill_histos or not os.path.exists(cacheFileName)
    optionsToPrint = ['inputDir', 'outputDir', 'region', 'tag', 'doFillHistograms']
    if verbose : print "options:\n"+'\n'.join(["%s : %s"%(o, eval(o)) for o in optionsToPrint])
    # collect inputs
    if verbose : print 'input files ',os.path.join(inputDir, templateInputFilename)
    tupleFilenames = glob.glob(os.path.join(inputDir, templateInputFilename))
    samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
    if not samples : samples = [guessSampleFromFilename(f) for f in tupleFilenames] # if the fast guess didn't work, try the slow one
    samplesPerGroup = collections.defaultdict(list)
    filenamesPerGroup = collections.defaultdict(list)
    for s, f in zip(samples, tupleFilenames) :
        samplesPerGroup[s.group].append(s)
        filenamesPerGroup[s.group].append(f)
    vars = ['pt','eta','d0sig','z0SinTheta','etCone','ptCone','etConeCorr','ptConeCorr']
    vars += ['relEtConeStd', 'relPtConeStd', 'relEtConeMod', 'relPtConeMod']
    groups = samplesPerGroup.keys()
    sources = leptonSources
    #fill histos
    if doFillHistograms :
        lepLabel = "(probe %s)"%lepton
        histosPerGroup = bookHistosPerGroup(vars, groups, lepLabel=lepLabel)
        histosPerSource = bookHistosPerSource(vars, sources, lepLabel=lepLabel)
        for group in groups:
            isData = isDataSample(group)
            filenames = filenamesPerGroup[group]
            histosThisGroup = histosPerGroup[group]
            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            print "%s : %d entries"%(group, chain.GetEntries())
            fillHistos(chain, histosThisGroup, histosPerSource, isData, lepton, group, verbose)
        writeHistos(cacheFileName, {'perGroup':histosPerGroup, 'perSource':histosPerSource}, verbose)
    # compute scale factors
    histosPerGroup = fetchHistos(cacheFileName, histoNames(vars, groups), verbose)
    histosPerSource = fetchHistos(cacheFileName, histoNames(vars, sources), verbose)
    plotStackedHistos(histosPerGroup,  outputDir+'/by_group',  region, colors=SampleUtils.colors, verbose=verbose)
    plotStackedHistos(histosPerSource, outputDir+'/by_source', region, colors=fakeu.colorsFillSources(), verbose=verbose)
    plotIsoComparison(histosPerSource, outputDir+'/',          region, lepton, verbose)

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

def bookHistosPerGroup(variables, groups, lepLabel='') :
    "book a dict of histograms with keys [sample][var][tight,loose]"
    def histo(variable, hname):
        h = None
        ptBinEdges = np.array([10.0, 20.0, 35.0, 60.0, 100.0]) # add a bin to see the denom effect above 60.0
        etaBinEdges = fakeu.etaBinEdges()
        if   v=='pt'          : h = r.TH1F(hname, ';p_{T}%s [GeV]; entries/bin'%lepLabel,      len(ptBinEdges)-1, ptBinEdges)
        elif v=='eta'         : h = r.TH1F(hname, ';#eta%s; entries/bin'%lepLabel,             len(etaBinEdges)-1, etaBinEdges)
        elif v=='d0sig'       : h = r.TH1F(hname, ';d_{0 sig}%s; entries/bin'%lepLabel,        50, -50.0, +50.0)
        elif v=='z0SinTheta'  : h = r.TH1F(hname, ';z_{0}sin#theta%s; entries/bin'%lepLabel,   50, -50.0, +50.0)
        elif v=='etCone'      : h = r.TH1F(hname, ';E_{T,cone}%s [GeV]; entries/bin'%lepLabel, 60, -10.0, +50.0)
        elif v=='ptCone'      : h = r.TH1F(hname, ';p_{T,cone}%s [GeV]; entries/bin'%lepLabel, 50,   0.0, +50.0)
        elif v=='etConeCorr'  : h = r.TH1F(hname, ';E_{T, cone, corr}%s; entries/bin'%lepLabel,60, -10.0, +50.0)
        elif v=='ptConeCorr'  : h = r.TH1F(hname, ';p_{T, cone, corr}%s; entries/bin'%lepLabel,50,   0.0, +50.0)
        elif v=='relEtConeStd': h = r.TH1F(hname, ';E_{T, cone, corr}/p_{T}%s; entries/bin'%lepLabel,        42, -0.02, +0.40)
        elif v=='relPtConeStd': h = r.TH1F(hname, ';p_{T, cone, corr}/p_{T}%s; entries/bin'%lepLabel,        42, -0.02, +0.40)
        elif v=='relEtConeMod': h = r.TH1F(hname, ';E_{T, cone, corr}/min(p_{T},60)%s; entries/bin'%lepLabel,42, -0.02, +0.40)
        elif v=='relPtConeMod': h = r.TH1F(hname, ';p_{T, cone, corr}/min(p_{T},60)%s; entries/bin'%lepLabel,42, -0.02, +0.40)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(g,
                  dict([(v,
                         {'loose'       : histo(variable=v, hname=histoNamePerSample(g, v, 'loose')),
                          'tight'       : histo(variable=v, hname=histoNamePerSample(g, v, 'tight')),
                          'tight_std'   : histo(variable=v, hname=histoNamePerSample(g, v, 'tight_std')),
                          'tight_tight' : histo(variable=v, hname=histoNamePerSample(g, v, 'tight_tight'))
                          })
                        for v in variables]))
                 for g in groups])

def bookHistosPerSource(variables, sources, lepLabel='') :
    "book a dict of histograms with keys [source][var][tight,loose]"
    return bookHistosPerGroup(variables=variables, groups=sources, lepLabel=lepLabel)


def extractName(dictOrHist):
    "input must be either a dict or something with 'GetName'"
    isDict = type(dictOrHist) is dict
    return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()

def histoNames(variables, samples) :
    return extractName(bookHistosPerGroup(variables, samples))

def fillHistos(chain, histosThisGroup, histosPerSource, isData, lepton, group, verbose=False):
    nLoose, nTight = 0, 0
    totWeightLoose, totWeightTight = 0.0, 0.0
    for event in chain :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        weight = weight
        probe = kin.addTlv(event.l1)
        met = kin.addTlv(event.met)
        isRightLep = probe.isEl if lepton=='el' else probe.isMu if lepton=='mu' else False
        if not isRightLep : continue
        # isTight = probe.isTight # this is the same as isTight_wh
        isTight       = fakeu.isTight_wh   (probe)
        isTight_std   = fakeu.isTight_std  (probe)
        isTight_tight = fakeu.isTight_tight(probe)
        if (isData and isTight != probe.isTight):
            print "isTight_wh %s isTight %s run %d event %d" % (isTight, probe.isTight, evtN, runN)
            print "(RunNumber==%d && EventNumber==%d)"%(runN, evtN)
            pt = probe.p4.Pt()
            relEtCone = probe.etConeCorr/min([pt, 60.0])
            relPtCone = probe.ptConeCorr/min([pt, 60.0])
            print "pt %.2f etCone %.2f relEt %.2f, ptCone %.2f relPt %.2f"%(pt, probe.etConeCorr, relEtCone, probe.ptConeCorr, relPtCone)
            fakeu.elecIsIsolatedVerbose(probe, fakeu.denominator_wh(probe))
        pt = probe.p4.Pt()
        eta = abs(probe.p4.Eta())
        d0Sig, z0Sin = probe.d0Signif, probe.z0SinTheta
        etCone, ptCone = probe.etCone, probe.ptCone
        etConeCorr, ptConeCorr = probe.etConeCorr, probe.ptConeCorr
        relDenomStd, relDenomMod = fakeu.denominator_std(probe), fakeu.denominator_wh(probe)
        relEtConeStd, relPtConeStd = etConeCorr*relDenomStd, ptConeCorr*relDenomStd
        relEtConeMod, relPtConeMod = etConeCorr*relDenomMod, ptConeCorr*relDenomMod
        def shiftWithinRange(v, maxV=0.40, epsilon=1.0e-3) : return maxV*(1.0-epsilon) if v>=maxV else v
        for vv in [relEtConeStd, relPtConeStd, relEtConeMod, relPtConeMod] : vv = shiftWithinRange(vv)
        relEtConeStd = shiftWithinRange(relEtConeStd)
        relPtConeStd = shiftWithinRange(relPtConeStd)
        relEtConeMod = shiftWithinRange(relEtConeMod)
        relPtConeMod = shiftWithinRange(relPtConeMod)
        selection = True
        source = None if isData else enum2source(probe)
        if selection:
            def incrementCounts(counts, weightedCounts):
                counts +=1
                weightedCounts += weight
            def fillPerGroup(tightOrLoose):
                histosThisGroup['pt'          ][tightOrLoose].Fill(pt, weight)
                histosThisGroup['eta'         ][tightOrLoose].Fill(eta, weight)
                histosThisGroup['d0sig'       ][tightOrLoose].Fill(d0Sig, weight)
                histosThisGroup['z0SinTheta'  ][tightOrLoose].Fill(z0Sin, weight)
                histosThisGroup['etCone'      ][tightOrLoose].Fill(etCone, weight)
                histosThisGroup['ptCone'      ][tightOrLoose].Fill(ptCone, weight)
                histosThisGroup['etConeCorr'  ][tightOrLoose].Fill(etConeCorr, weight)
                histosThisGroup['ptConeCorr'  ][tightOrLoose].Fill(ptConeCorr, weight)
                histosThisGroup['relEtConeStd'][tightOrLoose].Fill(relEtConeStd, weight)
                histosThisGroup['relPtConeStd'][tightOrLoose].Fill(relPtConeStd, weight)
                histosThisGroup['relEtConeMod'][tightOrLoose].Fill(relEtConeMod, weight)
                histosThisGroup['relPtConeMod'][tightOrLoose].Fill(relPtConeMod, weight)
                if source:
                    histosThisSource = histosPerSource[source]
                    histosThisSource['pt'          ][tightOrLoose].Fill(pt, weight)
                    histosThisSource['eta'         ][tightOrLoose].Fill(eta, weight)
                    histosThisSource['d0sig'       ][tightOrLoose].Fill(d0Sig, weight)
                    histosThisSource['z0SinTheta'  ][tightOrLoose].Fill(z0Sin, weight)
                    histosThisSource['etCone'      ][tightOrLoose].Fill(etCone, weight)
                    histosThisSource['ptCone'      ][tightOrLoose].Fill(ptCone, weight)
                    histosThisSource['etConeCorr'  ][tightOrLoose].Fill(etConeCorr, weight)
                    histosThisSource['ptConeCorr'  ][tightOrLoose].Fill(ptConeCorr, weight)
                    histosThisSource['relEtConeStd'][tightOrLoose].Fill(relEtConeStd, weight)
                    histosThisSource['relPtConeStd'][tightOrLoose].Fill(relPtConeStd, weight)
                    histosThisSource['relEtConeMod'][tightOrLoose].Fill(relEtConeMod, weight)
                    histosThisSource['relPtConeMod'][tightOrLoose].Fill(relPtConeMod, weight)

            incrementCounts(nLoose, totWeightLoose)
            fillPerGroup('loose')
            if isTight:
                incrementCounts(nTight, totWeightTight)
                fillPerGroup('tight')
            if isTight_std:
                fillPerGroup('tight_std')
            if isTight_tight:
                fillPerGroup('tight_tight')
    if verbose:
        counterNames = ['nLoose', 'nTight', 'totWeightLoose', 'totWeightTight']
        print ', '.join(["%s : %.1f"%(c, eval(c)) for c in counterNames])

def plotStackedHistos(histosPerGroup={}, outputDir='', region='', colors={}, verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    mkdirIfNeeded(outputDir)
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
        data = histosPerGroup['data'] if 'data' in histosPerGroup else None
        if data and data.GetEntries():
            data.SetMarkerStyle(r.kFullDotLarge)
            data.Draw('p same')
        yMin, yMax = getMinMax([h for h in [totBkg, data, err_band] if h])
        pm.SetMinimum(0.0)
        pm.SetMaximum(1.1*yMax)
        can.Update()
        topRightLabel(can, "#splitline{%s}{%s}"%(histoname, region), xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]+[data]
        can.Update()
        can.SaveAs(os.path.join(outputDir, histoname+'.png'))

def plotIsoComparison(histosPerSource={}, outputDir='', region='', lepton='', verbose=False):
    """
    plot a comparison of eff(T|L) for real and for fake leptons
    vs. pt, where the numerator is one of the tight definitions
    """
    var = 'pt'
    sources = histosPerSource.keys()
    lOrTOrTs = first(first(histosPerSource)).keys()
    histosPtPerSource = dict((s, dict((lt, histosPerSource[s][var][lt]) for lt in lOrTOrTs)) for s in sources)
    def buildTotFakeHistos():
        "add up all the non-real (fake) sources"
        notRealSources = [s for s in sources if s!='real']
        aSource = first(notRealSources)
        totFakeHistos = dict()
        for lt in ['loose', 'tight', 'tight_std', 'tight_tight']:
            template = histosPtPerSource[aSource][lt]
            h = template.Clone(template.GetName().replace(aSource, 'fake'))
            h.Reset()
            for s in sources : h.Add(histosPtPerSource[s][lt])
            totFakeHistos[lt] = h
        return totFakeHistos
    histosPtPerSource['fake'] = buildTotFakeHistos()
    effReal_wh    = rootUtils.buildRatioHistogram(histosPtPerSource['real']['tight'      ], histosPtPerSource['real']['loose'])
    effReal_std   = rootUtils.buildRatioHistogram(histosPtPerSource['real']['tight_std'  ], histosPtPerSource['real']['loose'])
    effReal_tight = rootUtils.buildRatioHistogram(histosPtPerSource['real']['tight_tight'], histosPtPerSource['real']['loose'])
    effFake_wh    = rootUtils.buildRatioHistogram(histosPtPerSource['fake']['tight'      ], histosPtPerSource['fake']['loose'])
    effFake_std   = rootUtils.buildRatioHistogram(histosPtPerSource['fake']['tight_std'  ], histosPtPerSource['fake']['loose'])
    effFake_tight = rootUtils.buildRatioHistogram(histosPtPerSource['fake']['tight_tight'], histosPtPerSource['fake']['loose'])
    frameName, frameTitle = region+'_'+lepton, "fake and real efficiencies for %s in %s"%(lepton, region)
    can = r.TCanvas('c_'+frameName, frameTitle, 800, 600)
    can.cd()
    pm = effReal_wh
    pm.SetMinimum(0.0)
    pm.SetMaximum(1.1)
    pm.GetYaxis().SetTitle("#epsilon(T|L)")
    colorReal, colorFake = r.kBlue, r.kRed
    markerWh, markerStd, markerTight = r.kMultiply, r.kCircle, r.kOpenSquare
    def setAttrs(h, mark, col):
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerStyle(mark)
    setAttrs(effReal_wh,    markerWh,    colorReal)
    setAttrs(effReal_std,   markerStd,   colorReal)
    setAttrs(effReal_tight, markerTight, colorReal)
    setAttrs(effFake_wh,    markerWh,    colorFake)
    setAttrs(effFake_std,   markerStd,   colorFake)
    setAttrs(effFake_tight, markerTight, colorFake)
    pm.SetStats(0)
    pm.Draw('axis')
    for h in [effReal_wh, effReal_std, effReal_tight, effFake_wh, effFake_std, effFake_tight]:
        h.Draw('same')
    leg = rightLegend(can)
    leg.SetBorderSize(0)
    leg.AddEntry(r.TObject(),   'Real', '')
    leg.AddEntry(effReal_std,   'std iso', 'lp')
    leg.AddEntry(effReal_tight, 'tight iso', 'lp')
    leg.AddEntry(effReal_wh,    'wh iso',  'lp')
    leg.AddEntry(r.TObject(),   'Fake', '')
    leg.AddEntry(effFake_std,   'std iso', 'lp')
    leg.AddEntry(effFake_tight, 'tight iso', 'lp')
    leg.AddEntry(effFake_wh,  '  wh iso',  'lp')
    leg.Draw()
    topRightLabel(can, "#splitline{%s}{%s}"%(lepton, region), xpos=0.125, align=13)
    can.RedrawAxis()
    can._histos = [effReal_wh, effReal_std, effFake_wh, effFake_std]
    can.Update()
    mkdirIfNeeded(outputDir)
    can.SaveAs(os.path.join(outputDir, frameTitle+'.png'))

if __name__=='__main__':
    main()
