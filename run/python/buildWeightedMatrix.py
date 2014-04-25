#!/bin/env python

# Combine the matrices for the matrix-method fake estimate as a
# weighted sum of the fake composition in each signal region.

# All the elements of the matrix, both p(T|R) and p(T|F), are in fact
# p_T-dependent histograms.
#
# For p(T|R) we simply count for each background the number of real
# leptons in each signal region.
# For p(T|F) we also subdivide the leptons in fake categories. That
# is, heavy-flavor for muons, heavy-flavor or conversion for
# electrons. The categorization as true/fake/hf/conv is based on the
# monte-carlo truth.
# After computing the fractions of leptons due to each background (and
# to each fake categories), a bin-by-bin weighted sum is performed,
# where the bins are now p_T-bins. That is, we perform a weighted sum
# of histograms.
#
# The procedure is described in sec. 6.2 of ATL-COM-PHYS-2012-1808.
# This python implementation is based on the c++ implementation by
# Matt (mrelich6@gmail.com), originally in FinalNewFake.cxx
#
# davide.gerbaudo@gmail.com
# October 2013

from math import sqrt
import operator
import optparse
import os
from rootUtils import importRoot, buildRatioHistogram, drawLegendWithDictKeys, getMinMax, getBinContents, getBinIndices, getBinning
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
from utils import (enumFromHeader
                   ,first
                   ,mkdirIfNeeded
                   ,json_write
                   ,rmIfExists
                   )
from fakeUtils import (samples
                       ,fakeProcesses
                       ,getInputFiles
                       ,buildRatio
                       )

import matplotlib as mpl
mpl.use('Agg') # render plots without X
import matplotlib.pyplot as plt
import numpy as np
import SampleUtils
from compute_fake_el_scale_factor import histoname_electron_sf_vs_eta

usage="""
Example usage:
%prog \\
 --tag ${TAG} \\
 --input_dir out/fakerate/merged/data_${TAG}.root \\
 --input-el-sf out/fake_el_scale_factor_${TAG}/hflf_el_scale_histos.root \\
 --input-el-sf out/fake_el_scale_factor_${TAG}/conv_el_scale_histos.root \\
 --output_file out/fakerate/merged/FinalFakeHist_${TAG}.root \\
 --output_plot out/fakerate/merged/FinalFakeHist_plots_${TAG} \\
 >& log/fakerate/FinalFakeHist_${TAG}.log
"""

# scale factors from determineFakeScaleFactor.py
# --- paste the lines below in buildWeightedMatrix.py ---
# Feb_12, 2014-02-12 18:20:20.650121
mu_qcdSF, mu_realSF = 0.86, 0.99590
el_convSF, el_qcdSF, el_realSF = 1.09, 0.63, 0.99633

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-f', '--input_fractions')
    parser.add_option('-o', '--output_file')
    parser.add_option('-p', '--output_plot')
    parser.add_option('-s', '--input-el-sf', default=[], action='append', help='electron bin-by-bin scale factors (from compute_fake_el_scale_factor)')
    parser.add_option('-z', '--zoom-in', help='vertical axis efficiency plots')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_dir', 'output_file', 'output_plot']
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag = opts.tag
    inputDirname  = opts.input_dir
    inputFracFname= opts.input_fractions
    inputSfFnames = opts.input_el_sf
    outputFname   = opts.output_file
    outputPlotDir = opts.output_plot
    zoomIn        = opts.zoom_in
    verbose       = opts.verbose
    if verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    allInputFiles = getInputFiles(inputDirname, tag, verbose) # includes allBkg, which is used only for sys
    assert all(f for f in allInputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in allInputFiles.iteritems()]))
    if inputSfFnames and any([not os.path.exists(f) for f in inputSfFnames]) : parser.error("invalid electron sf file(s) %s"%inputSfFnames)
    outputPlotDir = outputPlotDir+'/' if not outputPlotDir.endswith('/') else ''
    mkdirIfNeeded(outputPlotDir)
    outputFile = r.TFile.Open(outputFname, 'recreate')
    inputFiles = dict((k, v) for k, v in allInputFiles.iteritems() if k in fakeProcesses())
    inputFracFile = r.TFile.Open(inputFracFname) if inputFracFname else None
    if inputFracFname and not inputFracFile : parser.error("invalid fraction file %s"%inputFracFname)

    buildMuonRates    (inputFiles, outputFile, outputPlotDir, inputFracFile=inputFracFile, verbose=verbose, zoomIn=zoomIn)
    buildElectronRates(inputFiles, outputFile, outputPlotDir, inputFracFile=inputFracFile, inputElecSfFiles=inputSfFnames, verbose=verbose, zoomIn=zoomIn)
    buildSystematics  (allInputFiles['allBkg'], outputFile, verbose)
    outputFile.Close()
    if verbose : print "output saved to \n%s"%'\n'.join([outputFname, outputPlotDir])

def frac2str(frac) :
    flatFraction = type(first(frac)) is float
    return ('\n'.join([''.join("%12s"%s for s in fakeProcesses()),
                       ''.join("%12s"%("%.3f"%frac[s]) for s in fakeProcesses())])
            if flatFraction
            else '\n'.join([''.join("%12s"%s for s in fakeProcesses())]
                           +[''.join("%12s"%("%.3f"%frac[s].GetBinContent(b)) for s in fakeProcesses())
                             for b in getBinIndices(first(frac))]))

def selectionRegions() :
    print "hardcoded selectionRegions, should match what's in FakeRegions.h; fix DiLeptonMatrixMethod"
    return  ['CR_SRWHnoMlj', 'CR_SRWH1j', 'CR_SSInc1j']
#--tmp--    return ['CR_SSInc',
#--tmp--            'CR_SSInc1j',
#--tmp--            'CR_WHSS',
#--tmp--            'CR_CR8lpt',
#--tmp--            'CR_CR8ee',
#--tmp--            'CR_CR8mm',
#--tmp--            'CR_CR8mmMtww',
#--tmp--            'CR_CR8mmHt',
#--tmp--            'CR_CR9lpt',
#--tmp--            'CR_SsEwk',
#--tmp--            'CR_SsEwkLoose',
#--tmp--            'CR_SsEwkLea',
#--tmp--            'CR_WHZVfake1jee',
#--tmp--            'CR_WHZVfake2jee',
#--tmp--            'CR_WHZVfake1jem',
#--tmp--            'CR_WHZVfake2jem',
#--tmp--            'CR_WHfake1jem',
#--tmp--            'CR_WHfake2jem',
#--tmp--            'CR_WHZV1jmm',
#--tmp--            'CR_WHZV2jmm',
#--tmp--            'CR_WHfake1jmm',
#--tmp--            'CR_WHfake2jmm',
#--tmp--
#--tmp--            "CR_WHZVfake1j",
#--tmp--            "CR_WHZVfake2j",
#--tmp--            "CR_WHfake1j",
#--tmp--            "CR_WHfake2j",
#--tmp--            "CR_WHZV1j",
#--tmp--            "CR_WHZV2j",
#--tmp--
#--tmp--            "CR_SRWH1j",
#--tmp--            "CR_SRWH2j"
#--tmp--            ]

def isRegionToBePlotted(sr) :
    "regions for which we plot the weighted matrices"
    srs = ['CR_WHZVfake1j', 'CR_WHZVfake2j', 'CR_WHfake1j', 'CR_WHfake2j', 'CR_WHZV1j', 'CR_WHZV2j', 'CR_SRWH1j', 'CR_SRWH2j']
    srs = ['CR_SRWHnoMlj', 'CR_SRWH1j', 'CR_SSInc1j']
    return sr in srs

def extractionRegions() :
    return ['qcdMC', 'convMC', 'realMC']

def getRealEff(lepton='electron|muon', inputFile=None, scaleFactor=1.0) :
    histoName = lepton+'_realMC_all_l_pt_coarse'
    effHisto = buildRatio(inputFile, histoName)
    effHisto.Scale(scaleFactor)
    return effHisto
def buildRatioAndScaleIt(histoPrefix='', inputFile=None, scaleFactor=1.0, verbose=False) :
    ratioHisto = buildRatio(inputFile, histoPrefix)
    def lf2s(l) : return ', '.join(["%.3f"%e for e in l])
    if verbose: print ratioHisto.GetName()," before scaling: ",lf2s(getBinContents(ratioHisto))
    if   type(scaleFactor)==float : ratioHisto.Scale(scaleFactor)
    elif type(scaleFactor)==r.TH1F : ratioHisto.Multiply(scaleFactor)
    elif type(scaleFactor)==r.TH2F : ratioHisto.Multiply(scaleFactor)
    else : raise TypeError("unknown SF type %s, %s"%(type(scaleFactor), str(scaleFactor)))
    if verbose: print ratioHisto.GetName()," after scaling: ",lf2s(getBinContents(ratioHisto))
    return ratioHisto
def buildPercentages(inputFiles, histoName, binLabel) :
    "build a dictionary (process, fraction of counts) for a given bin"
    histos = dict((p, f.Get(histoName)) for p, f in inputFiles.iteritems())
    assert all(h for h in histos.values()),"missing histo '%s'\n%s"%(histoName, '\n'.join("%s : %s"%(k, str(v)) for k, v in histos.iteritems()))
    bin = first(histos).GetXaxis().FindBin(binLabel)
    counts = dict((p, h.GetBinContent(bin)) for p, h in histos.iteritems())
    norm = sum(counts.values())
    if not norm : print "buildPercentages: warning, all empty histograms for %s[%s]"%(histoName, binLabel)
    counts = dict((p, c/norm if norm else 0.0)for p,c in counts.iteritems())
    return counts
def buildPercentagesTwice(inputFiles, histoName, binLabelA, binLabelB) :
    """
    build two dictionaries (process, fraction of counts) for two given bins.
    Note that we cannot use buildPercentages because we want to
    normalize A and B (i.e. qcd and conv) together.
    Maybe refactor these two buildPercentages* functions?
    """
    histos = dict((p, f.Get(histoName)) for p, f in inputFiles.iteritems())
    assert all(h for h in histos.values()),"missing histo '%s'\n%s"%(histoName, '\n'.join("%s : %s"%(k, str(v)) for k, v in histos.iteritems()))
    binA = first(histos).GetXaxis().FindBin(binLabelA)
    binB = first(histos).GetXaxis().FindBin(binLabelB)
    countsA = dict((p, h.GetBinContent(binA)) for p, h in histos.iteritems())
    countsB = dict((p, h.GetBinContent(binB)) for p, h in histos.iteritems())
    norm = sum(countsA.values() + countsB.values())
    if not norm : print "buildPercentages: warning, all empty histograms for %s[%s,%s]"%(histoName, binLabelA, binLabelB)
    countsA = dict((p, c/norm if norm else 0.0)for p,c in countsA.iteritems())
    countsB = dict((p, c/norm if norm else 0.0)for p,c in countsB.iteritems())
    return countsA, countsB
def binWeightedSum(histos={}, weights={}, bin=1) :
    assert not set(histos)-set(weights), "different keys: histos[%s], weights[%s]"%(str(histos.keys()), str(weights.keys()))
    cews = [(h.GetBinContent(bin), h.GetBinError(bin), w)
            for h, w in [(histos[k], weights[k]) for k in histos.keys()]]
    tot  = sum(c*w     for c, e, w in cews)
    err2 = sum(e*e*w*w for c, e, w in cews)
    print "weighted bin %d (w=%.1f) : %.3f = %s"%(bin, sum(w for c, e, w in cews), tot,
                                                  '+'.join("%.2f*%.2f"%(w,c) for c, e, w in cews))
    return tot, err2
def buildWeightedHisto(histos={}, fractions={}, histoName='', histoTitle='') :
    "was getFinalRate"
    hout = first(histos).Clone(histoName if histoName else 'final_rate') # should pick a better default
    hout.SetTitle(histoTitle)
    hout.Reset()
    flatFraction = type(first(fractions)) is float
    if flatFraction :
        print "averaging flat ",histoName
        print 'keys -> ',histos.keys()
        for b in getBinIndices(hout) :
            tot, err2 = binWeightedSum(histos, fractions, b)
            hout.SetBinContent(b, tot)
            hout.SetBinError(b, sqrt(err2))
    else :
        bH, bF = getBinning(first(histos)), getBinning(first(fractions))
        assert bH == bF,"different binning %s: %s, %s: %s"%(first(histos).GetName(), bH, first(fractions).GetName(), bF)
        weightedHistos = dict((p, h.Clone(h.GetName()+'_weighted_for_'+histoName)) for p,h in histos.iteritems()) # preserve originals
        print "averaging 2d ",histoName
        for b in getBinIndices(hout):
            print "bin %d (w=%.1f):  %.3f = %s"%(b,
                                                 sum(fractions[p].GetBinContent(b) for p in fractions.keys()),
                                                 sum(fractions[p].GetBinContent(b)*weightedHistos[p].GetBinContent(b) for p in fractions.keys()),
                                                 '+'.join("%.2f*%.2f"%(fractions[p].GetBinContent(b), weightedHistos[p].GetBinContent(b))
                                                          for p in fractions.keys()))
        for p,h in weightedHistos.iteritems() :
            h.Multiply(fractions[p])
            hout.Add(h)
    return hout
def buildWeightedHistoTwice(histosA={}, fractionsA={}, histosB={}, fractionsB={},
                            histoName='', histoTitle='') :
    "was getFinalRate"
    assert not set(histosA)-set(histosB),"different keys A[%s], B[%s]"%(str(histosA.keys()), str(histosB.keys()))
    assert type(first(fractionsA)) is type(first(fractionsB)),"etherogenous fractions A %s, B %s"%(type(first(fractionsA)), type(first(fractionsB)))
    hout = first(histosA).Clone(histoName if histoName else 'final_rate') # should pick a better default
    hout.SetTitle(histoTitle)
    hout.Reset()
    flatFraction = type(first(fractionsA)) is float
    if flatFraction :
        print "averaging flat ",histoName
        print 'keys -> ',histosA.keys()
        for b in getBinIndices(hout) :
            totA, errA2 = binWeightedSum(histosA, fractionsA, b)
            totB, errB2 = binWeightedSum(histosB, fractionsB, b)
            hout.SetBinContent(b, totA + totB)
            hout.SetBinError(b, sqrt(errA2 + errB2))
    else :
        bHA, bFA = getBinning(first(histosA)), getBinning(first(fractionsA))
        bHB, bFB = getBinning(first(histosB)), getBinning(first(fractionsB))
        assert bHA == bFA,"different binning %s: %s, %s: %s"%(first(histosA).GetName(), bHA, first(fractionsA).GetName(), bFA)
        assert bHB == bFB,"different binning %s: %s, %s: %s"%(first(histosB).GetName(), bHB, first(fractionsB).GetName(), bFB)
        assert bHA == bHB,"different binning %s: %s, %s: %s"%(first(histosA).GetName(), bHA, first(histosB).GetName(), bHB)
        weightedHistosA = dict((p, h.Clone(h.GetName()+'_weighted_for_'+histoName)) for p,h in histosA.iteritems()) # preserve originals
        weightedHistosB = dict((p, h.Clone(h.GetName()+'_weighted_for_'+histoName)) for p,h in histosB.iteritems())


        print "averaging 2d ",histoName
        for b in getBinIndices(hout):
            print "bin %d (w=%.1f):  %.3f = %s"%(b,
                                                 sum(fractionsA[p].GetBinContent(b)+fractionsB[p].GetBinContent(b) for p in fractionsA.keys()),
                                                 sum(fractionsA[p].GetBinContent(b)*weightedHistosA[p].GetBinContent(b)
                                                     +fractionsB[p].GetBinContent(b)*weightedHistosB[p].GetBinContent(b)
                                                     for p in fractionsA.keys()),
                                                 '+'.join(["%.2f*%.2f"%(fractionsA[p].GetBinContent(b),
                                                                        weightedHistosA[p].GetBinContent(b))
                                                           for p in fractionsA.keys()]+
                                                          ["%.2f*%.2f"%(fractionsB[p].GetBinContent(b),
                                                                        weightedHistosB[p].GetBinContent(b))
                                                           for p in fractionsB.keys()]))


        for p in weightedHistosA.keys() :
            hA, fA = weightedHistosA[p], fractionsA[p]
            hB, fB = weightedHistosB[p], fractionsB[p]
            hA.Multiply(fA)
            hB.Multiply(fB)
            hout.Add(hA)
            hout.Add(hB)
    return hout
def buildMuonRates(inputFiles, outputfile, outplotdir, inputFracFile=None, verbose=False, zoomIn=False) :
    """
    For each selection region, build the real eff and fake rate
    histo as a weighted sum of the corresponding fractions.
    """
    processes = fakeProcesses()
    brsit, iF, v = buildRatioAndScaleIt, inputFiles, verbose
    print "buildMuonRates: values to be fixed: ",' '.join(["%s: %s"%(_, eval(_)) for _ in ['mu_qcdSF', 'mu_realSF']])
    eff_qcd  = dict((p, brsit('muon_qcdMC_all_l_pt_coarse',  iF[p], mu_qcdSF, v))  for p in processes)
    eff_real = dict((p, brsit('muon_realMC_all_l_pt_coarse', iF[p], mu_realSF, v)) for p in processes)
    eff2d_qcd  = dict((p, brsit('muon_qcdMC_all_l_pt_eta',  iF[p], mu_qcdSF, v))  for p in processes)
    eff2d_real = dict((p, brsit('muon_realMC_all_l_pt_eta', iF[p], mu_realSF, v)) for p in processes)
    lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#varepsilon(T|L)'
    plot1dEfficiencies(eff_qcd,  'eff_mu_qcd',  outplotdir, lT+' qcd fake #mu'+';'+lX+';'+lY, zoomIn=zoomIn)
    plot1dEfficiencies(eff_real, 'eff_mu_real', outplotdir, lT+' real #mu'    +';'+lX+';'+lY, zoomIn=zoomIn)
    lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#eta'
    plot2dEfficiencies(eff2d_qcd,  'eff2d_mu_qcd', outplotdir, lT+' qcd fake #mu'+';'+lX+';'+lY)
    plot2dEfficiencies(eff2d_real, 'eff2d_mu_real', outplotdir, lT+' real #mu'   +';'+lX+';'+lY)
    mu_frac = dict()
    for sr in selectionRegions() :
        fC, bP = fetchCompositions, buildPercentages
        isf = inputFracFile
        hnTemplateQcd = '%(proc)s_muon_'+sr+'_all_flavor_pt_den_qcd'
        hnTemplateReal = '%(proc)s_muon_'+sr+'_all_flavor_pt_den_real'
        frac_qcd  = fC(isf, hnTemplateQcd,  processes) if isf else bP(inputFiles, 'muon_'+sr+'_all_flavor_den', 'qcd')
        frac_real = fC(isf, hnTemplateReal, processes) if isf else bP(inputFiles, 'muon_'+sr+'_all_flavor_den', 'real')
        if verbose : print "mu : sr ",sr,"\n frac_qcd  : ",frac2str(frac_qcd )
        if verbose : print "mu : sr ",sr,"\n frac_real : ",frac2str(frac_real)
        fake1d = buildWeightedHisto(eff_qcd,  frac_qcd, 'mu_fake_rate_'+sr, 'Muon fake rate '+sr)
        real1d = buildWeightedHisto(eff_real, frac_real, 'mu_real_eff_'+sr, 'Muon real eff ' +sr)

        c2dC = compose2Dcompositions
        hnTemplateQcd = '%(proc)s_muon_'+sr+'_all_flavor_pt_%(etabin)s_den_qcd'
        hnTemplateReal = '%(proc)s_muon_'+sr+'_all_flavor_pt_%(etabin)s_den_real'
        frac_qcd2d  = c2dC(isf, hnTemplateQcd,  processes) if isf else frac_qcd
        frac_real2d = c2dC(isf, hnTemplateReal, processes) if isf else frac_real
        fake2d = buildWeightedHisto(eff2d_qcd,  frac_qcd2d, 'mu_fake_rate2d_'+sr, 'Muon fake rate #eta vs. p_{T}'+sr)
        real2d = buildWeightedHisto(eff2d_real, frac_real2d, 'mu_real_eff2d_'+sr, 'Muon real eff  #eta vs. p_{T}'+sr)

        outputfile.cd()
        fake1d.Write()
        real1d.Write()
        fake2d.Write()
        real2d.Write()
        mu_frac[sr] = {'qcd' : frac_qcd, 'real' : frac_real}
        if isRegionToBePlotted(sr) :
            print "plotting eff2d_mu_fake:"
            fake2d.Print('all')
            print "plotting eff2d_mu_real"
            real2d.Print('all')
            lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#eta'
            plot2dEfficiencies({sr : fake2d}, 'eff2d_mu_fake', outplotdir, lT+' fake #mu'+';'+lX+';'+lY)
            plot2dEfficiencies({sr : real2d}, 'eff2d_mu_real', outplotdir, lT+' real #mu'+';'+lX+';'+lY)
    #json_write(mu_frac, outplotdir+/outFracFilename)
    doPlotFractions = not inputFracFile
    if doPlotFractions : plotFractions(mu_frac, outplotdir, 'frac_mu')
def fetchSfHistos(inputElecSfFiles=[], histoname='', verbose=False):
    # inputElecSfFile should be two files, one for conv and one for bbcc
    fileNames = inputElecSfFiles if type(inputElecSfFiles)==list else inputElecSfFiles.split()
    assert len(fileNames)==2,"fetchSfHistos expects two files (hflf+conv), got %s"%str(inputElecSfFile)
    if verbose : print "retrieving scale factors from %s"%inputElecSfFiles
    fname_hflf = first(filter(lambda _ : 'hflf' in _, fileNames))
    fname_conv = first(filter(lambda _ : 'conv' in _, fileNames))
    file_hflf = r.TFile.Open(fname_hflf)
    file_conv = r.TFile.Open(fname_conv)
    hname = histoname_electron_sf_vs_eta()
    histos = {'hflf' : composeEtaHistosAs2dPtEta(input1Dhisto=file_hflf.Get(hname), outhistoname=hname+'_hflf'),
              'conv' : composeEtaHistosAs2dPtEta(input1Dhisto=file_conv.Get(hname), outhistoname=hname+'_conv')
              }
    for f in [file_hflf, file_conv] : f.Close()
    return histos
def buildElectronRates(inputFiles, outputfile, outplotdir, inputFracFile=None, inputElecSfFiles=[], verbose=False, zoomIn=False) :
    """
    For each selection region, build the real eff and fake rate
    histo as a weighted sum of the corresponding fractions.
    Note that the fake has two components (conversion and qcd).
    """
    processes = fakeProcesses()
    brsit, iF, v = buildRatioAndScaleIt, inputFiles, verbose
    # if we have the sf vs eta, use it in the 2d parametrization
    el_convSF_vs_eta = el_convSF if not inputElecSfFiles else fetchSfHistos(inputElecSfFiles, histoname_electron_sf_vs_eta, verbose)['conv']
    el_qcdSF_vs_eta  = el_qcdSF if not inputElecSfFiles else fetchSfHistos(inputElecSfFiles, histoname_electron_sf_vs_eta, verbose)['hflf']
    if inputElecSfFiles : el_convSF_vs_eta.Print("all")
    if inputElecSfFiles : el_qcdSF_vs_eta.Print("all")

    eff_conv = dict((p, brsit('elec_convMC_all_l_pt_coarse', iF[p], el_convSF, v)) for p in processes)
    eff_qcd  = dict((p, brsit('elec_qcdMC_all_l_pt_coarse',  iF[p], el_qcdSF, v))  for p in processes)
    eff_real = dict((p, brsit('elec_realMC_all_l_pt_coarse', iF[p], el_realSF, v)) for p in processes)
    eff2d_conv = dict((p, brsit('elec_convMC_all_l_pt_eta', iF[p], el_convSF_vs_eta, v)) for p in processes)
    eff2d_qcd  = dict((p, brsit('elec_qcdMC_all_l_pt_eta',  iF[p], el_qcdSF_vs_eta, v))  for p in processes)
    eff2d_real = dict((p, brsit('elec_realMC_all_l_pt_eta', iF[p], el_realSF, v)) for p in processes)

    lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#varepsilon(T|L)'
    plot1dEfficiencies(eff_conv, 'eff_el_conv', outplotdir, lT+' conv fake el'+';'+lX+';'+lY, zoomIn=zoomIn)
    plot1dEfficiencies(eff_qcd,  'eff_el_qcd',  outplotdir, lT+' qcd fake el' +';'+lX+';'+lY, zoomIn=zoomIn)
    plot1dEfficiencies(eff_real, 'eff_el_real', outplotdir, lT+' real el'     +';'+lX+';'+lY, zoomIn=zoomIn)
    lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#eta'
    plot2dEfficiencies(eff2d_conv, 'eff2d_el_conv', outplotdir, lT+' conv fake el'+';'+lX+';'+lY)
    plot2dEfficiencies(eff2d_qcd,  'eff2d_el_qcd',  outplotdir, lT+' qcd fake el' +';'+lX+';'+lY)
    plot2dEfficiencies(eff2d_real, 'eff2d_el_real', outplotdir, lT+' real el'     +';'+lX+';'+lY)

    el_frac = dict()
    for sr in selectionRegions() :
        fC, bPt = fetchCompositions, buildPercentagesTwice
        isf = inputFracFile
        hnTemplateQcd  = '%(proc)s_elec_'+sr+'_all_flavor_pt_den_qcd'
        hnTemplateConv = '%(proc)s_elec_'+sr+'_all_flavor_pt_den_conv'
        hnTemplateReal = '%(proc)s_elec_'+sr+'_all_flavor_pt_den_real'
        frac_conv, frac_qcd= (fC(isf, hnTemplateConv,  processes), fC(isf, hnTemplateQcd,  processes)) if isf else bPt(inputFiles, 'elec_'+sr+'_all_flavor_den', 'conv', 'qcd')
        frac_real = fC(isf, hnTemplateReal, processes) if isf else buildPercentages(inputFiles, 'elec_'+sr+'_all_flavor_den', 'real')
        if verbose : print "el : sr ",sr,"\n frac_conv : ",frac2str(frac_conv)
        if verbose : print "el : sr ",sr,"\n frac_qcd  : ",frac2str(frac_qcd )
        if verbose : print "el : sr ",sr,"\n frac_real : ",frac2str(frac_real)
        real1d = buildWeightedHisto     (eff_real, frac_real,                     'el_real_eff_'+sr, 'Electron real eff '+sr)
        fake1d = buildWeightedHistoTwice(eff_conv, frac_conv, eff_qcd,  frac_qcd, 'el_fake_rate_'+sr, 'Electron fake rate '+sr)

        c2dC = compose2Dcompositions
        hnTemplateQcd  = '%(proc)s_elec_'+sr+'_all_flavor_pt_%(etabin)s_den_qcd'
        hnTemplateConv = '%(proc)s_elec_'+sr+'_all_flavor_pt_%(etabin)s_den_conv'
        hnTemplateReal = '%(proc)s_elec_'+sr+'_all_flavor_pt_%(etabin)s_den_real'
        frac_qcd2d  = c2dC(isf, hnTemplateQcd,  processes) if isf else frac_qcd
        frac_conv2d = c2dC(isf, hnTemplateConv, processes) if isf else frac_conv
        frac_real2d = c2dC(isf, hnTemplateReal, processes) if isf else frac_real
        real2d = buildWeightedHisto     (eff2d_real, frac_real2d,                     'el_real_eff2d_'+sr, 'Electron real eff  #eta vs. p_{T}'+sr)
        fake2d = buildWeightedHistoTwice(eff2d_conv, frac_conv2d, eff2d_qcd,  frac_qcd2d, 'el_fake_rate2d_'+sr, 'Electron fake rate  #eta vs. p_{T}'+sr)

        outputfile.cd()
        fake1d.Write()
        real1d.Write()
        fake2d.Write()
        real2d.Write()
        el_frac[sr] = {'conv' : frac_conv, 'qcd' : frac_qcd, 'real' : frac_real}
        if isRegionToBePlotted(sr) :
            lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#eta'
            plot2dEfficiencies({sr : fake2d}, 'eff2d_el_fake', outplotdir, lT+' fake e'+';'+lX+';'+lY)
            plot2dEfficiencies({sr : real2d}, 'eff2d_el_real', outplotdir, lT+' real e'+';'+lX+';'+lY)
    #json_write(el_frac, outFracFilename)
    doPlotFractions = not inputFracFile
    if doPlotFractions : plotFractions(el_frac, outplotdir, 'frac_el')
def buildEtaSyst(inputFileTotMc, inputHistoBaseName='(elec|muon)_qcdMC_all', outputHistoName='', verbose=False) :
    """
    Take the eta distribution and normalize it to the average fake
    rate (taken from one bin rate); use the differences from 1 as the
    fractional uncertainty.
    """
    rate = buildRatio(inputFileTotMc, inputHistoBaseName+'_l_eta_coarse').Clone(outputHistoName)
    norm = buildRatio(inputFileTotMc, inputHistoBaseName+'_onebin').GetBinContent(1)
    rate.Scale(1.0/norm if norm else 1.0)
    bins = range(1, 1+rate.GetNbinsX())
    for b in bins : rate.AddBinContent(b, -1.0) # DG there must be a better way to do this
    scaleUpForward, fwFact, maxCentralEta = True, 2.0, 1.5
    if scaleUpForward :
        for b in bins :
            bCon, bCen = rate.GetBinContent(b), rate.GetBinCenter(b)
            rate.SetBinContent(b, bCon*(fwFact if abs(bCen)>maxCentralEta else 1.0))
    if inputHistoBaseName.startswith('mu') : rate.Reset() # mu consistent with 0.
    if verbose : print "eta syst ",inputHistoBaseName," : ",["%.2f"%rate.GetBinContent(b) for b in range(1, 1+rate.GetNbinsX())]
    return rate
def buildSystematics(inputFileTotMc, outputfile, verbose=False) :
    "Hardcoded values from FinalNewFake.h; might not be used at all...ask Matt"
    print "build syst might be droppped...check this with Matt"
    print "rename *_down to *_do"
    el_real_up = r.TParameter('double')('el_real_up', 0.01)
    el_real_dn = r.TParameter('double')('el_real_down', 0.02)
    mu_real_up = r.TParameter('double')('mu_real_up', 0.00)
    mu_real_dn = r.TParameter('double')('mu_real_down', 0.02)
    el_HFLFerr = r.TParameter('double')('el_HFLFerr', 0.05)
    mu_HFLFerr = r.TParameter('double')('mu_HFLFerr', 0.00)
    el_datamc  = r.TParameter('double')('el_datamc',  0.20) #datamc are effectively the sf error.
    mu_datamc  = r.TParameter('double')('mu_datamc',  0.05) #Right now taking the Pt variation into account
    el_region  = r.TParameter('double')('el_region',  0.05)
    mu_region  = r.TParameter('double')('mu_region',  0.10)
    el_eta     = buildEtaSyst(inputFileTotMc, 'elec_qcdMC_all', 'el_eta_sys', verbose)
    mu_eta     = buildEtaSyst(inputFileTotMc, 'muon_qcdMC_all', 'mu_eta_sys', verbose)
    allSys = [el_real_up, el_real_dn, mu_real_up, mu_real_dn,
              el_HFLFerr, mu_HFLFerr, el_datamc , mu_datamc,
              el_region, mu_region, el_eta, mu_eta ]
    outputfile.cd()
    for o in  allSys : o.Write()
def plotFractions(fractDict={}, outplotdir='./', prefix='') :
    """
    input : fractDict[sr][lep_type][sample] = float
    """
    outplotdir = outplotdir if outplotdir.endswith('/') else outplotdir+'/'
    #def isRegionToBePlotted(r) : return r in selectionRegions()+extractionRegions()
    regions  = sorted(filter(isRegionToBePlotted, fractDict.keys()))
    leptypes = sorted(first(fractDict).keys())
    samples  = sorted(first(first(fractDict)).keys())
    ind = np.arange(len(regions))
    width = 0.5
    colors = dict(zip(samples, ['b','g','r','c','m','y']))
    for lt in leptypes :
        fracPerSample = dict((s, np.array([fractDict[r][lt][s] for r in regions])) for s in samples)
        below = np.zeros(len(regions))
        plots = []
        fig, ax = plt.subplots()
        for s, frac in fracPerSample.iteritems() :
            plots.append(plt.bar(ind, frac, width, color=colors[s], bottom=below))
            below = below + frac
        plt.ylabel('fractions')
        plt.title(prefix+' '+lt+' compositions')
        plt.xticks(ind+width/2., regions)
        plt.ylim((0.0, 1.0))
        plt.grid(True)
        plt.yticks(np.arange(0.0, 1.0, 0.2))
        labels = {'heavyflavor' : 'bb/cc', 'diboson' : 'VV', 'ttbar':'tt'}
        labels = [labels[s] if s in labels else s for s in samples]
        leg = plt.legend([p[0] for p in plots], labels, bbox_to_anchor=(1.135, 1.05))
        leg.get_frame().set_alpha(0.5)
        fig.autofmt_xdate(bottom=0.25, rotation=90, ha='center')
        fig.savefig(outplotdir+prefix+'_'+lt+'.png')
        fig.savefig(outplotdir+prefix+'_'+lt+'.eps')

def plot1dEfficiencies(effs={}, canvasName='', outputDir='./', frameTitle='title;p_{T} [GeV]; efficiency', zoomIn=False) :
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    padMaster = None
    colors, markers = SampleUtils.colors, SampleUtils.markers
    for s,h in effs.iteritems() :
        h.SetLineColor(colors[s] if s in colors else r.kBlack)
        h.SetMarkerColor(h.GetLineColor())
        h.SetMarkerStyle(markers[s] if s in markers else r.kFullCircle)
        drawOpt = 'ep same' if padMaster else 'ep'
        h.Draw(drawOpt)
        if not padMaster : padMaster = h
    minY, maxY = getMinMax(effs.values()) if zoomIn else (0.0, 1.0)
    padMaster.GetYaxis().SetRangeUser(min([0.0, minY]), 1.1*maxY)
    padMaster.SetMinimum(0.0)
    padMaster.SetMaximum(1.1*maxY)
    padMaster.SetMaximum(0.25)
    padMaster.SetTitle(frameTitle)
    padMaster.SetStats(False)
    drawLegendWithDictKeys(can, effs)
    can.Update()
    for ext in ['png','eps'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)
def plot2dEfficiencies(effs={}, canvasName='', outputDir='./', frameTitle='efficiency; #eta; p_{T} [GeV]', zoomIn=False) :
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    origTextFormat = r.gStyle.GetPaintTextFormat()
    r.gStyle.SetPaintTextFormat('.2f')
    for s,h in effs.iteritems() :
        can.Clear()
        # todo minZ, maxZ = getMinMax(effs.values()) if zoomIn else (0.0, 1.0)
        minZ, maxZ = (0.0, 1.0)
        h.SetMarkerSize(1.5*h.GetMarkerSize())
        h.Draw('colz')
        h.Draw('text e same')
        h.GetZaxis().SetRangeUser(min([0.0, minZ]), maxZ)
        def dropCrPrefix(sr) : return sr.replace('CR_', '')
        h.SetTitle(dropCrPrefix(s)+' : '+frameTitle)
        h.SetStats(False)
        can.Update()
        for ext in ['png','eps'] :
            outFilename = outputDir+'/'+canvasName+'_'+s+'.'+ext
            rmIfExists(outFilename)
            can.SaveAs(outFilename)
    r.gStyle.SetPaintTextFormat(origTextFormat)
def fetchCompositions(inputSfFile=None, templateHistoName="%(proc)s", processes=[]) :
    histos = dict((p, inputSfFile.Get(templateHistoName%{'proc':p})) for p in processes)
    missing = dict([(k, v) for k, v in histos.iteritems() if not v])
    assert not missing,"missing compositions: (template %s) from %s:\n%s"%(templateHistoName, inputSfFile.GetName(), str(missing))
    return histos
def compose2Dcompositions(inputSfFile=None, templateHistoName="%(proc)s_%(etabin)s", processes=[]) :
    "take two 1D fractions histograms for one eta slice each, and compose them in a 2D fractions histogram"
    etaBins = ['etaC', 'etaF']
    histos1d = dict((e, dict((p, inputSfFile.Get(templateHistoName%{'proc':p, 'etabin':e})) for p in processes)) for e in etaBins)
    assert all(v for ve in histos1d.values() for v in ve.values()),"missing compositions: %s"%histos1d
    h1dC, h1dF = first(histos1d['etaC']), first(histos1d['etaF'])
    nX, xMin, xMax = h1dC.GetNbinsX(), h1dC.GetXaxis().GetBinLowEdge(1), h1dC.GetXaxis().GetBinUpEdge(h1dC.GetNbinsX())
    histos2d = dict((p, r.TH2F(histos1d['etaC'][p].GetName().replace('_etaC_','_vs_eta_'), '', nX, xMin, xMax, 2, 0.0, 2.0)) for p in histos1d['etaC'].keys())
    for p in processes :
        for iEta, eta in zip(range(1, 1+len(etaBins)), etaBins) :
            hEta, hEtaPt = histos1d[eta][p], histos2d[p]
            for iPt in range(1, 1+nX) :
                hEtaPt.SetBinContent(iPt, iEta, hEta.GetBinContent(iPt))
                hEtaPt.SetBinError  (iPt, iEta, hEta.GetBinError  (iPt))
    return histos2d

def composeEtaHistosAs2dPtEta(input1Dhisto=None, outhistoname='') :
    "take the 1D scale factor histogram (vs eta), and build a 2D histo that has (pt,eta) on (x,y); see MeasureFakeRate2::initHistos"
    ptBinEdges = np.array([10.0, 20.0, 35.0, 100.0]) # see FakeBinnings.h -> coarseFakePtbins
    etaBinEdges = np.array([0.0, 1.37, 2.50])
    h = r.TH2F(outhistoname, '', len(ptBinEdges)-1, ptBinEdges, len(etaBinEdges)-1, etaBinEdges)
    h.SetDirectory(0)
    h.Sumw2()
    for iX in range(1, 1+len(ptBinEdges)):
        for iY in range(1, 1+len(etaBinEdges)):
            h.SetBinContent(iX, iY, input1Dhisto.GetBinContent(iY))
            h.SetBinError  (iX, iY, input1Dhisto.GetBinError  (iY))
    return h

if __name__=='__main__' :
    main()
