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
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from rootUtils import buildRatioHistogram
from utils import (enumFromHeader
                   ,first
                   ,json_write
                   )
import matplotlib as mpl
mpl.use('Agg') # render plots without X
import matplotlib.pyplot as plt
import numpy as np

usage="""
Example usage:
%prog \\
 --tag ${TAG} \\
 --input_dir out/fakerate/merged/data_${TAG}.root \\
 --output_file out/fakerate/merged/FinalFakeHist_${TAG}.root \\
 --output_plot out/fakerate/merged/FinalFakeHist_plots_${TAG} \\
 >& log/fakerate/FinalFakeHist_${TAG}.log
"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-o', '--output_file')
    parser.add_option('-p', '--output_plot')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_dir', 'output_file', 'output_plot']
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag = opts.tag
    inputDirname  = opts.input_dir
    outputFname   = opts.output_file
    outputPlotDir = opts.output_plot
    verbose       = opts.verbose
    if verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    allInputFiles = getInputFiles(inputDirname, tag, verbose) # includes allBkg, which is used only for sys
    assert all(f for f in allInputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in allInputFiles.iteritems()]))
    outputFile = r.TFile.Open(outputFname, 'recreate')
    inputFiles = dict((k, v) for k, v in allInputFiles.iteritems() if k in fakeProcesses())

    buildMuonRates    (inputFiles, outputFile, outputPlotDir, verbose)
    buildElectronRates(inputFiles, outputFile, outputPlotDir, verbose)
    buildSystematics  (allInputFiles['allBkg'], outputFile)
    outputFile.Close()
    if verbose : print "output saved to \n%s"%'\n'.join([outputFname, outFractionsMu, outFractionsEl])

def samples() : return ['allBkg', 'ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def fakeProcesses() : return ['ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def frac2str(frac) :
    return '\n'.join([''.join("%12s"%s for s in fakeProcesses()),
                      ''.join("%12s"%("%.3f"%frac[s]) for s in fakeProcesses())])
def selectionRegions() :
    header = os.path.dirname(__file__)+'/../../SusyTest0/SusyAnaDefsMatt.h'
    enum = enumFromHeader(header, 'SignalRegion')
    enum = [x[0] for x in sorted(enum.iteritems(), key=operator.itemgetter(1))] # sort by value
    enum = enum[:-1] # the last one is just the enum size (SR_N), not an actual value
    fix = {'CRSSInc':'CR_SSInc',
           'SR_WHSS':'CR_WHSS',
           'CR8lpt'    :'CR_CR8lpt'    ,
           'CR8ee'     :'CR_CR8ee'     ,
           'CR8mm'     :'CR_CR8mm'     ,
           'CR8mmMtww' :'CR_CR8mmMtww' ,
           'CR8mmHt'   :'CR_CR8mmHt'   ,
           } # some enums have string repr different from their literal repr
    print "selectionRegions : fix ugly mapping when conflicting enums have been removed"
    enum = [fix[e] if e in fix else e for e in enum]
    return enum
def getInputFiles(inputDirname, tag, verbose=False) :
    inDir = inputDirname
    tag = tag if tag.startswith('_') else '_'+tag
    files = dict(zip(samples(), [r.TFile.Open(inDir+'/'+s+tag+'.root') for s in samples()]))
    if verbose : print "getInputFiles('%s'):\n\t%s"%(inputDirname, '\n\t'.join("%s : %s"%(k, f.GetName()) for k, f in files.iteritems()))
    return files
def buildRatio(inputFile=None, histoBaseName='') :
    num, den = inputFile.Get(histoBaseName+'_num'), inputFile.Get(histoBaseName+'_den')
    return buildRatioHistogram(num, den, histoBaseName +'_rat')
def getRealEff(lepton='electron|muon', inputFile=None, scaleFactor=1.0) :
    histoName = lepton+'_realMC_all_l_pt_coarse'
    effHisto = buildRatio(inputFile, histoName)
    effHisto.Scale(scaleFactor)
    return effHisto
def buildRatioAndScaleIt(histoPrefix='', inputFile=None, scaleFactor=1.0) :
    ratioHisto = buildRatio(inputFile, histoPrefix)
    ratioHisto.Scale(scaleFactor)
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
    return tot, err2
def buildWeightedHisto(histos={}, fractions={}, histoName='', histoTitle='') :
    "was getFinalRate"
    hout = first(histos).Clone(histoName if histoName else 'final_rate') # should pick a better default
    hout.SetTitle(histoTitle)
    hout.Reset()
    for b in range(1, 1+hout.GetNbinsX()) :
        tot, err2 = binWeightedSum(histos, fractions, b)
        hout.SetBinContent(b, tot)
        hout.SetBinError(b, sqrt(err2))
    return hout
def buildWeightedHistoTwice(histosA={}, fractionsA={}, histosB={}, fractionsB={},
                            histoName='', histoTitle='') :
    "was getFinalRate"
    assert not set(histosA)-set(histosB),"different keys A[%s], B[%s]"%(str(histosA.keys()), str(histosB.keys()))
    hout = first(histosA).Clone(histoName if histoName else 'final_rate') # should pick a better default
    hout.SetTitle(histoTitle)
    hout.Reset()
    for b in range(1, 1+hout.GetNbinsX()) :
        totA, errA2 = binWeightedSum(histosA, fractionsA, b)
        totB, errB2 = binWeightedSum(histosB, fractionsB, b)
        hout.SetBinContent(b, totA + totB)
        hout.SetBinError(b, sqrt(errA2 + errB2))
    return hout
def buildMuonRates(inputFiles, outputfile, outplotdir, verbose=False) :
    """
    For each selection region, build the real eff and fake rate
    histo as a weighted sum of the corresponding fractions.
    """
    processes = fakeProcesses()
    brsit, iF = buildRatioAndScaleIt, inputFiles
    mu_qcdSF, mu_realSF = 0.79059, 0.99719
    print "buildMuonRates: values to be fixed: ",' '.join(["%s: %s"%(v, eval(v)) for v in ['mu_qcdSF', 'mu_realSF']])
    eff_qcd  = dict((p, brsit('muon_qcdMC_all_l_pt_coarse',  iF[p], mu_qcdSF))  for p in processes)
    eff_real = dict((p, brsit('muon_realMC_all_l_pt_coarse', iF[p], mu_realSF)) for p in processes)
    mu_frac = dict()
    for sr in selectionRegions() :
        frac_qcd  = buildPercentages(inputFiles, 'muon_'+sr+'_all_flavor_den', 'qcd')
        frac_real = buildPercentages(inputFiles, 'muon_'+sr+'_all_flavor_den', 'real')
        if verbose : print "mu : sr ",sr,"\n frac_qcd  : ",frac2str(frac_qcd )
        if verbose : print "mu : sr ",sr,"\n frac_real : ",frac2str(frac_real)
        fake = buildWeightedHisto(eff_qcd,  frac_qcd, 'mu_fake_rate_'+sr, 'Muon fake rate '+sr)
        real = buildWeightedHisto(eff_real, frac_real, 'mu_real_eff_'+sr, 'Muon real eff '+sr)
        outputfile.cd()
        fake.Write()
        real.Write()
        mu_frac[sr] = {'qcd' : frac_qcd, 'real' : frac_real}
    #json_write(mu_frac, outplotdir+/outFracFilename)
    plotFractions(mu_frac, outplotdir, 'mu')
def buildElectronRates(inputFiles, outputfile, outplotdir, verbose=False) :
    """
    For each selection region, build the real eff and fake rate
    histo as a weighted sum of the corresponding fractions.
    Note that the fake has two components (conversion and qcd).
    """
    processes = fakeProcesses()
    brsit, iF = buildRatioAndScaleIt, inputFiles
    el_convSF, el_qcdSF, el_realSF = 1.24359, 0.73345, 0.99729
    print "buildElectronRates: values to be fixed: ",' '.join(["%s: %s"%(v, eval(v)) for v in ['el_qcdSF', 'el_convSF', 'el_realSF']])
    eff_conv = dict((p, brsit('elec_convMC_all_l_pt_coarse', iF[p], el_convSF)) for p in processes)
    eff_qcd  = dict((p, brsit('elec_qcdMC_all_l_pt_coarse',  iF[p], el_qcdSF))  for p in processes)
    eff_real = dict((p, brsit('elec_realMC_all_l_pt_coarse', iF[p], el_realSF)) for p in processes)
    el_frac = dict()
    for sr in selectionRegions() :
        frac_conv, frac_qcd= buildPercentagesTwice(inputFiles, 'elec_'+sr+'_all_flavor_den',
                                                   'conv', 'qcd')
        frac_real = buildPercentages(inputFiles, 'elec_'+sr+'_all_flavor_den', 'real')
        if verbose : print "el : sr ",sr,"\n frac_conv : ",frac2str(frac_conv)
        if verbose : print "el : sr ",sr,"\n frac_qcd  : ",frac2str(frac_qcd )
        if verbose : print "el : sr ",sr,"\n frac_real : ",frac2str(frac_real)
        real = buildWeightedHisto     (eff_real, frac_real,                     'el_real_eff_'+sr, 'Electron real eff '+sr)
        fake = buildWeightedHistoTwice(eff_conv, frac_conv, eff_qcd,  frac_qcd, 'el_fake_rate_'+sr, 'Electron fake rate '+sr)
        outputfile.cd()
        fake.Write()
        real.Write()
        el_frac[sr] = {'conv' : frac_conv, 'qcd' : frac_qcd, 'real' : frac_real}
    #json_write(el_frac, outFracFilename)
    plotFractions(el_frac, outplotdir, 'el')
def buildEtaSyst(inputFileTotMc, inputHistoBaseName='(elec|muon)_qcdMC_all', outputHistoName='') :
    """
    Take the eta distribution and normalize it to the average fake
    rate (taken from one bin rate); use the differences from 1 as the
    fractional uncertainty.
    """
    rate = buildRatio(inputFileTotMc, inputHistoBaseName+'_l_eta_coarse').Clone(outputHistoName)
    norm = buildRatio(inputFileTotMc, inputHistoBaseName+'_onebin').GetBinContent(1)
    rate.Scale(1.0/norm if norm else 1.0)
    for b in range(1, 1+rate.GetNbinsX()) : rate.AddBinContent(b, -1.0) # DG there must be a better way to do this
    if inputHistoBaseName.startswith('mu') : rate.Reset() # mu consistent with 0.
    return rate
def buildSystematics(inputFileTotMc, outputfile) :
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
    el_eta     = buildEtaSyst(inputFileTotMc, 'elec_qcdMC_all', 'el_eta_sys')
    mu_eta     = buildEtaSyst(inputFileTotMc, 'muon_qcdMC_all', 'mu_eta_sys')
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
    def isInterestingRegion(r) : return any(k in r for k in ['CR8', 'WHSS', 'SSInc'])
    regions  = [r for r in selectionRegions() if isInterestingRegion(r)]
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
        plt.legend([p[0] for p in plots], samples, bbox_to_anchor=(1.2, 1.05))
        fig.autofmt_xdate(bottom=0.25, rotation=90, ha='center')
        plt.savefig(outplotdir+prefix+'_'+lt+'.png')

if __name__=='__main__' :
    main()
