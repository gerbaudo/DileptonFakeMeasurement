#!/bin/env python

# script to build the fake weighted matrices from the fake ntuples

# This script puts together the various pieces determined from the
# fake nutples. Inputs:
# - compute_eff_from_ntuple.py
# - compute_fake_compositions.py
# - compute_fake_el_scale_factor.py
#
# output: weighted average fake matrices (for now electron only)
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
from compute_fake_el_scale_factor import histoname_sf_vs_eta

usage="""
Example usage:
%prog \\
 --verbose  \\
 --output-dir ./out/fakerate/el_sf_${TAG}
 >& log/fakerate/el_sf_${TAG}.log

 TODO
"""
# old flat values (still use them for 1D parametrization?)
el_convSF, el_qcdSF = 1.09, 1.3
mu_qcdSF = 1.3

def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-c', '--comp-histos', help='output from compute_fake_compositions.py')
    parser.add_option('-e', '--eff-histos', default=[], action='append', help='output files from compute_eff_from_ntuple.py')
    parser.add_option('-r', '--region', help='where we have the compositions, and want the fake matrix, e.g. ssinc1j, emu')
    parser.add_option('-s', '--scale-factors', default=[], action='append', help='bin-by-bin data/mc from compute_fake_scale_factor')
    parser.add_option('-o', '--output-dir', default='./out/fake_weighted_average', help='dir for plots')
    parser.add_option('-l', '--lepton', default='el', help='either el or mu')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('--also-anygroup', action='store_true', help='also build matrix without compositions,'                      ' to evaluate the systematic uncertainty on the composition')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    (options, args) = parser.parse_args()
    inputDir  = options.input_dir
    compFname = options.comp_histos
    effFnames = options.eff_histos
    region    = options.region
    sfFnames  = options.scale_factors
    outputDir = options.output_dir
    lepton    = options.lepton
    verbose   = options.verbose
    if lepton not in ['el', 'mu'] : parser.error("invalid lepton '%s'"%lepton)
    if region not in ['emu', 'ssinc', 'ssinc1j','razor0j'] : parser.error("invalid region '%s'"%region)
    if not compFname or not os.path.exists(compFname) : parser.error("invalid composition file '%s'"%compFname)
    if not effFnames or not all(os.path.exists(f) for f in effFnames) : parser.error("invalid efficiency file '%s'"%str(effFnames))
    if not sfFnames  or not all(os.path.exists(f) for f in sfFnames) :
        # parser.error("invalid electron sf file(s) %s"%str(sfFnames))
        print "missing sf files, using flat scale factors"# do not crash, fall back on flat scale factors
    optionsToPrint = ['inputDir', 'outputDir']
    if verbose :
        print "working from %s"%os.getcwd()
        print "being called as : %s"%' '.join(os.sys.argv)
        print "options parsed:\n"+'\n'.join(["%s : %s"%(o, eval(o)) for o in optionsToPrint])
    mkdirIfNeeded(outputDir)
    # collect inputs
    regions = [region, ]
    groups=['diboson', 'heavyflavor', 'ttbar', 'wjets', 'zjets']
    if lepton=='el' : groups = filter(lambda _ : _!='heavyflavor', groups) # note this must be in sync with compute_fake_compositions; TODO: implement a way to get the groups from the available histos
    compositions = fetchCompositionHistos(compFname, lepton, groups, regions, verbose) # [var][group][reg][orig], note here orig=[conv,heavy,light]
    # pprint.pprint(compositions)
    efficiencies = fetchEffienciesHistos(effFnames, lepton, groups, verbose) # [var][group][orig], note here orig=[conv,heavy,light,qcd]
    # pprint.pprint(efficiencies)
    scale_factor_histos = fetchSfHistos(sfFnames, lepton, verbose)
    convSF_vs_eta = scale_factor_histos['conv'] if lepton=='el' else None
    qcdSF_vs_eta  = scale_factor_histos['hflf'] if 'hflf' in scale_factor_histos else None
    # if verbose:
    #     print "convSF: "+("vs. eta {0}".format(getBinContents(convSF_vs_eta)) if convSF_vs_eta else el_convSF)
    #     print "qcdSF: "+("vs. eta {0}".format(getBinContents(qcdSF_vs_eta)) if qcdSF_vs_eta else el_qcdSF)
    def scale_factor_to_str(sf):
        if 'vs_eta' in sf: return '['+', '.join("%.3f" % _ for _ in sf['vs_eta'])+']'
        else : return "%.3f" % sf['flat']
    if lepton=='el':
        scaleFactors = {'conv' : {'flat' : el_convSF, 'vs_eta' : convSF_vs_eta},
                        'heavy' : {'flat' : el_qcdSF, 'vs_eta' : qcdSF_vs_eta} }
        if verbose : print_scale_factor_dict(scaleFactors)
        scaleFakeEfficiencies(efficiencies, scaleFactors)
    elif lepton=='mu':
        scaleFactors = {'heavy' : {'flat' : mu_qcdSF, 'vs_eta' : qcdSF_vs_eta} }
        if verbose : print_scale_factor_dict(scaleFactors)
        scaleFakeEfficiencies(efficiencies, scaleFactors)

    # for now compute the weighted avg only for 'ssinc1j'
    avgEfficiencies = dict()
    for reg in first(first(compositions)).keys():
        avgEfficiencies[reg] = dict()
        for var in ['pt', 'pt_eta']:
            is1D = var=='pt'
            lT = "%s #varepsilon(T|L) fake %s"%(reg, lepton)
            lX = 'p_{T} [GeV]'
            lY = '#varepsilon(T|L)' if is1D else '#eta'
            hname = "%(lep)s_fake_%(var)s_%(reg)s"%{'lep':lepton, 'var':var, 'reg':reg}
            htitle = lT+';'+lX+';'+lY
            groups  = first(compositions).keys()
            origins = first(first(first(compositions))).keys()
            if verbose : print 'origins :',origins,'\n' + 'groups :',groups
            histosEff  = dict((group+'_'+orig, efficiencies[var][group]     [orig]) for group in groups for orig in origins)
            histosComp = dict((group+'_'+orig, compositions[var][group][reg][orig]) for group in groups for orig in origins)
            avgEff =  weightedAverage(histosEff, histosComp, hname, htitle, verbose)
            avgEfficiencies[reg][var] = avgEff
            if is1D:
                fakeu.plot1dEfficiencies({reg : avgEff}, 'eff1d_'+lepton+'_fake_'+reg, outputDir, htitle, zoomIn=True)
            else:
                fakeu.plot2dEfficiencies({reg : avgEff}, 'eff2d_'+lepton+'_fake_'+reg, outputDir, htitle, zoomIn=True)
    writeHistos(os.path.join(outputDir,'fake_matrices_'+lepton+'.root'), avgEfficiencies, verbose)
    if options.also_anygroup:# test with the group-independent efficiencies
	print 'fetchCompositionHistos ',compFname
	compositions = fetchCompositionHistos(compFname, lepton, ['anygroup'], verbose)
	pprint.pprint(compositions)
	print 'fetchEffienciesHistos ',effFnames
	efficiencies = fetchEffienciesHistos(effFnames, lepton, ['anygroup'], verbose)
	avgEfficiencies = dict()
	for reg in first(first(compositions)).keys():
	    avgEfficiencies[reg] = dict()
	    for var in ['pt', 'pt_eta']:
	        is1D = var=='pt'
	        lT = "%s #varepsilon(T|L) fake %s"%(reg, lepton)
	        lX = 'p_{T} [GeV]'
	        lY = '#varepsilon(T|L)' if is1D else '#eta'
	        hname = "%(lep)s_fake_%(var)s_%(reg)s"%{'lep':lepton, 'var':var, 'reg':reg}
	        htitle = lT+';'+lX+';'+lY
	        groups  = first(compositions).keys()
	        origins = first(first(first(compositions))).keys()
	        if verbose : print 'origins :',origins,'\n' + 'groups :',groups
	        histosEff  = dict((group+'_'+orig, efficiencies[var][group]     [orig])
                                  for group in groups for orig in origins)
	        histosComp = dict((group+'_'+orig, compositions[var][group][reg][orig])
                                  for group in groups for orig in origins)
	        avgEff =  weightedAverage(histosEff, histosComp, hname, htitle, verbose)
	        avgEfficiencies[reg][var] = avgEff
	        if is1D:
	            fakeu.plot1dEfficiencies({reg : avgEff}, 'eff1d_'+lepton+'_fake_'+reg+'_anygroup',
                                             outputDir, htitle, zoomIn=True)
	        else:
	            fakeu.plot2dEfficiencies({reg : avgEff}, 'eff2d_'+lepton+'_fake_'+reg+'_anygroup',
                                             outputDir, htitle, zoomIn=True)
	writeHistos(os.path.join(outputDir,'fake_matrices_'+lepton+'_anygroup.root'),
                    avgEfficiencies, verbose)

#___________________________________________________

leptonTypes = fakeu.leptonTypes()
allLeptonSources = fakeu.allLeptonSources()
leptonSources = fakeu.leptonSources()
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source
fetchSfHistos = fakeu.fetchSfHistos

def histoname_electron_sf_vs_eta() : return 'sf_el_vs_eta'
def histoname_electron_sf_vs_pt() : return 'sf_el_vs_pt'


def fetchCompositionHistos(filename, lepton, groups=[], regions=[], verbose=False):
    template = "h_%(var)s_%(group)s_%(region)s_%(origin)s_loose"
    origins = ['heavy', 'light']
    origins += ['conv'] if lepton=='el' else []
    #origins += ['real'] # real later, for now just worry about the hf/lf splitting
    variables = ['pt', 'pt_eta']
    histonames = dict((v,
                       dict((g,
                             dict((reg,
                                   dict((o,
                                         template%{'var':v, 'group':g, 'region':reg, 'origin':o})
                                        for o in origins))
                                  for reg in regions))
                            for g in groups))
                      for v in variables)
    return fetchHistos(filename, histonames, verbose)
def fetchEffienciesHistos(filenames, lepton,
                          groups=[],
                          verbose=False):
    template = "h_%(var)s_%(group)s_%(origin)s_tight_over_loose_%(region)s"
    regions = ['conv', 'hflf',] if lepton=='el' else ['hflf'] # these are the regions where we extract the eff, see MeasureFakeRate2
    fnameQcd  = first(filter(lambda _ : 'mcqcd' in _, filenames))
    fnameConv = first(filter(lambda _ : 'mcconv' in _, filenames))
    assert all(f and os.path.exists(f) for f in [fnameQcd, fnameConv] if f),"inputs: %s, qcd: %s, conv %s"%(str(filenames), fnameQcd, fnameConv)
    originsConv, originsQcd = ['conv'], ['heavy', 'light', 'qcd']
    vars = ['pt', 'pt_eta']
    def histonamesDict(origins=[], region=''):
        return dict((v,
                     dict((g, dict((o, template%{'var':v, 'group':g, 'region':region, 'origin':o}) for o in origins))
                          for g in groups))
                    for v in vars)
    histosConv = fetchHistos(fnameConv, histonamesDict(originsConv, 'conv'), verbose)
    histosQcd  = fetchHistos(fnameQcd,  histonamesDict(originsQcd,  'hflf'), verbose)
    def mergeConvQcd():
        return dict((v, dict((g, dict((o, histosConv[v][g][o] if o in originsConv else histosQcd[v][g][o])
                                      for o in originsQcd+originsConv)) for g in groups)) for v in vars)
    return mergeConvQcd() if lepton=='el' else histosQcd

def writeHistos(outputFileName='', histosPerSamplePerSource={}, verbose=False):
    rootUtils.writeObjectsToFile(outputFileName, histosPerSamplePerSource, verbose)

def fetchHistos(fileName='', histoNames={}, verbose=False):
    return rootUtils.fetchObjectsFromFile(fileName, histoNames, verbose)

def scaleFakeEfficiencies(efficiencies={}, scaleFactors={}):
    for var in efficiencies.keys():
        for group in first(efficiencies).keys():
            def is1d(h): return 'TH1' in h.ClassName()
            def is2d(h): return 'TH2' in h.ClassName()
            def scaleHisto(h, sf):
                if is1d(h) or not 'vs_eta' in sf or not sf['vs_eta'] :   h.Scale(sf['flat'])
                elif is2d(h) : h.Multiply(sf['vs_eta'])
            scaleHisto(efficiencies[var][group]['heavy'], scaleFactors['heavy'])
            scaleHisto(efficiencies[var][group]['qcd'], scaleFactors['heavy'])
            if 'conv' not in efficiencies[var][group] : continue # mu does not have a conv SF, we're done
            scaleHisto(efficiencies[var][group]['conv'], scaleFactors['conv'])

def weightedAverage(histosEff={}, histosWeight={}, histoName='', histoTitle='', verbose=False):
    getBinIndices, getBinContents, getBinning = rootUtils.getBinIndices, rootUtils.getBinContents, rootUtils.getBinning
    assert sorted(histosEff.keys())==sorted(histosWeight.keys()),"effs and weights must have the same keys:\n\teffs %s\n\tweights %s"%(str(histosEff.keys()), str(histosWeight.keys()))
    hout = first(histosEff).Clone(histoName if histoName else 'weighted_average_eff')
    hout.SetTitle(histoTitle)
    hout.Reset()
    def validateWeights(h):
        bins, values = getBinIndices(h), getBinContents(h)
        allValid = all(v>=0.0 and v<=1.0 for v in values)
        if not allValid:
            if verbose :
                print "warning '%s' weights not in [0,1] : [%s]"%(h.GetName(), ', '.join(("%.3f"%v for v in values)))
                print "setting them to 0.0 or 1.0"
            for b, v in zip(bins, values) : h.SetBinContent(b, 0.0 if v<0.0 else 1.0 if v>1.0 else v)
        return allValid
    [validateWeights(hw) for hw in histosWeight.values()]
    bins, binning = getBinIndices(first(histosEff)), getBinning(first(histosEff))
    for h in histosWeight.values()+histosEff.values():
        if getBinning(h)!=binning : print "warning: %s has binning %s, expecting %s"%(h.GetName(), str(getBinning(h)), str(binning))
    groups = sorted(histosEff.keys())
    epsilon = 1.0e-3
    binWeightNormalizations = [sum(histosWeight[g].GetBinContent(b) for g in groups) for b in bins]
    weightsAreNormalized = all(abs(1.0-norm)<epsilon for norm in binWeightNormalizations)
    if not weightsAreNormalized : print "warning, compositions are not normalized : [%s]"%', '.join(("%.3f"%v for v in binWeightNormalizations))
    print '-- ',histoName,'-- '
    for g in groups:
        bws, bcs = getBinContents(histosWeight[g]), getBinContents(histosEff[g])
        print "adding %18s : %s"%(g, ' : '.join("%.4f*%.4f"%(bw, bc) for bw, bc in zip(bws, bcs)))
        histoEff, histoWeight = histosEff[g], histosWeight[g]
        histoEff.Multiply(histoWeight)
        hout.Add(histoEff)
    print "tot weight   : %s"%' '.join(("%.4f"%v for v in (sum(histosWeight[g].GetBinContent(b) for g in groups) for b in bins)))
    print "weighted avg : %s"%' '.join(("%.4f"%v for v in getBinContents(hout)))
    return hout

def print_scale_factor_dict(sf_dict):
    "dict is something like {conv,heavy : {flat,vs_eta : values}}"
    print "scale factors:"
    for k,v in sf_dict.iteritems():
        values = getBinContents(v['vs_eta']) if 'vs_eta' in v and v['vs_eta'] else [v['flat'],]
        formatted_values = ', '.join("%.3f"%_ for _ in values)
        print "%s : %s"%(k, formatted_values)

if __name__=='__main__':
    main()
