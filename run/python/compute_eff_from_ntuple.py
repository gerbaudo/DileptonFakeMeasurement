#!/bin/env python

# script to compute the fake and real efficiencies from the MeasureFakeRate2 ntuples (as opposed to histos)

# davide.gerbaudo@gmail.com
# May 2014

import array
import collections
import glob
import math
import numpy as np
import optparse
import os
import pprint
import time
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
import kin
import fakeUtils as fakeu
import plotParametrizedFractions

usage="""
Example usage:
%prog \\
 --verbose  \\
 --tag ${TAG} \\
 --lepton el \\
 --mode conv \\
 --input-dir ./out/fakerate \\
 --output-dir ./out/fakerate/efficiencies_${TAG} \\
 >& log/fakerate/eff_el_conv_${TAG}.log
"""
def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='./out/fakerate/efficiencies')
    parser.add_option('-l', '--lepton', default='el', help='either el or mu')
    parser.add_option('-m', '--mode', help='real, conv, hflf')
    parser.add_option('-t', '--tag', help='tag used to select the input files (e.g. Apr_04)')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('-T', '--tight-def', help='on-the-fly tight def, one of defs in fakeUtils.py: fakeu.lepIsTight_std, etc.')
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
    validModesEl = ['real', 'hflf'] + ['conv']
    validModesMu = ['real', 'hflf']
    if mode not in (validModesEl if lepton=='el' else validModesMu) : parser.error("invalid mode %s"%mode)
    tupleStem, treeName = {'conv' : ('mcconv_tuple', 'ConversionExtractionRegion'),
                           'hflf' : ('mcqcd_tuple', 'HfLfExtractionRegion'),
                           'real' : ('mcreal_tuple', 'RealExtractionRegion')
                           }[mode]
    templateInputFilename = "*_%(stem)s_%(tag)s.root" % {'tag':tag, 'stem':tupleStem}
    templateOutputFilename =  "%(stem)s_%(l)s_eff.root" % {'stem':tupleStem.replace('tuple','histos'), 'l':lepton}
    outputFileName = os.path.join(outputDir, templateOutputFilename)
    cacheFileName = outputFileName.replace('.root', '_'+mode+'_cache.root')
    doFillHistograms = options.fill_histos or not os.path.exists(cacheFileName)
    onthefly_tight_def = eval(options.tight_def) if options.tight_def else None # eval will take care of aborting on typos
    optionsToPrint = ['inputDir', 'outputDir', 'mode', 'tag', 'doFillHistograms', 'cacheFileName', 'onthefly_tight_def']
    if verbose :
        print "working from %s"%os.getcwd()
        print "being called as : %s"%' '.join(os.sys.argv)
        print "options parsed:\n"+'\n'.join(["%s : %s"%(o, eval(o)) for o in optionsToPrint])
    # collect inputs
    print 'input filenames: ',os.path.join(inputDir, templateInputFilename)

    tupleFilenames = glob.glob(os.path.join(inputDir, templateInputFilename))
    samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
    samplesPerGroup = collections.defaultdict(list)
    filenamesPerGroup = collections.defaultdict(list)
    mkdirIfNeeded(outputDir)
    for s, f in zip(samples, tupleFilenames) :
        samplesPerGroup[s.group].append(s)
        filenamesPerGroup[s.group].append(f)
    vars = ['pt', 'pt_eta']
    groups = [g for g in samplesPerGroup.keys() if not isDataSample(g) and not g=='higgs']
    if lepton=='el' : groups = [g for g in groups if g!='heavyflavor']

    sourcesThisMode = {'real' : ['real'], # use same convention as in FakeLeptonSources.h
                       'conv' : ['conv'],
                       'hflf' : ['heavy', 'light', 'qcd']
                       }[mode]
    #fill histos
    if doFillHistograms :
        start_time = time.clock()
        num_processed_entries = 0
        histosPerGroupPerSource = bookHistosPerSamplePerSource(vars, groups, sourcesThisMode, mode=mode)
        for group in groups:
            filenames = filenamesPerGroup[group]
            histosThisGroupPerSource = dict((v, histosPerGroupPerSource[v][group]) for v in histosPerGroupPerSource.keys())
            histosAnyGroupPerSource  = dict((v, histosPerGroupPerSource[v]['anygroup']) for v in histosPerGroupPerSource.keys())

            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            if verbose: print "%s : %d entries"%(group, chain.GetEntries())
            num_processed_entries += fillHistos(chain, histosThisGroupPerSource, histosAnyGroupPerSource,
                                                lepton, mode,
                                                onthefly_tight_def=onthefly_tight_def, verbose=verbose)
        writeHistos(cacheFileName, histosPerGroupPerSource, verbose)
        end_time = time.clock()
        delta_time = end_time - start_time
        if verbose:
            print ("processed {0:d} entries ".format(num_processed_entries)
                   +"in "+("{0:d} min ".format(int(delta_time/60)) if delta_time>60 else
                           "{0:.1f} s ".format(delta_time))
                   +"({0:.1f} kHz)".format(num_processed_entries/delta_time))
    # compute efficiencies
    histosPerGroupPerSource = fetchHistos(cacheFileName, histoNamesPerSamplePerSource(vars, groups, sourcesThisMode, mode), verbose)
    effs = computeEfficiencies(histosPerGroupPerSource) # still [var][gr][source][l/t]
    for s in sourcesThisMode:
        for v in vars:
            groups = first(effs).keys()
            varIs1D, varIs2D = v=='pt', v=='pt_eta'
            effsThisSourceThisVar = dict((g, effs[v][g][s]) for g in groups)
            densThisSourceThisVar = dict((g, histosPerGroupPerSource[v][g][s]['loose']) for g in groups if g!='anygroup')
            numsThisSourceThisVar = dict((g, histosPerGroupPerSource[v][g][s]['tight']) for g in groups if g!='anygroup')
            if varIs1D:
                cname = 'eff_'+lepton+'_'+s
                lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#varepsilon(T|L)'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                zoomIn = True
                fakeu.plot1dEfficiencies(effsThisSourceThisVar, cname, outputDir, title, zoomIn)
                cname = 'stack_loose_'+lepton+'_'+s
                lT, lY = 'loose '+lepton+', denominator to #varepsilon(T|L)', '#varepsilon(T|L)'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                plotParametrizedFractions.plotStackedHistos(densThisSourceThisVar, cname, outputDir, title)
                cname = 'stack_tight_'+lepton+'_'+s
                lT, lY = 'tight '+lepton+', numerator to #varepsilon(T|L)', '#varepsilon(T|L)'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                plotParametrizedFractions.plotStackedHistos(numsThisSourceThisVar, cname, outputDir, title)

            elif varIs2D:
                cname = 'eff_'+lepton+'_'+s
                lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#eta'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                fakeu.plot2dEfficiencies(effsThisSourceThisVar, cname, outputDir, title, zoomIn=zoomIn)
    writeHistos(outputFileName, effs, verbose)
    if verbose : print "saved scale factors to %s" % outputFileName

#___________________________________________________

leptonTypes = ['loose', 'tight']
leptonSources = []
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def fillHistos(chain, histosPerSource, histosPerSourceAnygroup, lepton, mode, onthefly_tight_def=None, verbose=False):
    "fill the histograms, returns the number of events processed"
    class Counters: # scope trick (otherwise unavailable within nested func
        nLoose, nTight = 0, 0
        totWeightLoose, totWeightTight = 0.0, 0.0
        def str(self):
            counterNames = ['nLoose', 'nTight', 'totWeightLoose', 'totWeightTight']
            return ', '.join(["%s : %.1f"%(c, getattr(self, c)) for c in counterNames])
    counters = Counters()
    addTlv = kin.addTlv
    hspsAnygroup = histosPerSourceAnygroup
    if verbose : print "about to loop on {0} entries".format(chain.GetEntries())
    num_processed_entries = 0
    for iEvent, event in enumerate(chain) :
        num_processed_entries += 1
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        l0, l1 = event.l0, event.l1
        def fillHistosBySource(lep):
            lep = addTlv(lep)
            isTight = onthefly_tight_def(lep) if onthefly_tight_def else lep.isTight
            source = enum2source(lep)
            isRightLep = lep.isEl if lepton=='el' else lep.isMu
            isRightSource = (mode=='real' and source=='real' or
                             mode=='conv' and source=='conv' or
                             mode=='hflf' and source in ['heavy', 'light'])
            def fill(tightOrLoose):
                if tightOrLoose=='loose':
                    counters.nLoose +=1
                    counters.totWeightLoose += weight
                if tightOrLoose=='tight':
                    counters.nTight +=1
                    counters.totWeightTight += weight
                histosPerSource['pt'    ][source][tightOrLoose].Fill(pt,      weight)
                histosPerSource['pt_eta'][source][tightOrLoose].Fill(pt, eta, weight)
                hspsAnygroup   ['pt'    ][source][tightOrLoose].Fill(pt,      weight)
                hspsAnygroup   ['pt_eta'][source][tightOrLoose].Fill(pt, eta, weight)
                if mode=='hflf': # qcd is heavy+light
                    histosPerSource['pt'    ]['qcd'][tightOrLoose].Fill(pt,      weight)
                    histosPerSource['pt_eta']['qcd'][tightOrLoose].Fill(pt, eta, weight)
                    hspsAnygroup   ['pt'    ]['qcd'][tightOrLoose].Fill(pt,      weight)
                    hspsAnygroup   ['pt_eta']['qcd'][tightOrLoose].Fill(pt, eta, weight)
            if isRightLep and isRightSource :
                lep = addTlv(lep)
                pt, eta = lep.p4.Pt(), abs(lep.p4.Eta())
                fill('loose')
                if isTight : fill('tight')
        fillHistosBySource(l0)
        fillHistosBySource(l1)
    if verbose : print counters.str()
    return num_processed_entries

def histoNamePerSamplePerSource(var, sample, leptonSource, tightOrLoose, mode):
    return 'h_'+var+'_'+sample+'_'+leptonSource+'_'+tightOrLoose+'_'+mode
def bookHistosPerSamplePerSource(variables, samples, sources, mode=''):
    "book a dict of histograms with keys [var][sample][lepton_source][tight, loose]"
    print 'booking histos'
    def histo(var, hname):
        h = None
        ptBinEdges = fakeu.ptBinEdges()
        etaBinEdges = fakeu.etaBinEdges()
        if   var=='pt'     : h = r.TH1F(hname, ';p_{T} [GeV]',       len(ptBinEdges)-1,  ptBinEdges)
        elif var=='pt_eta' : h = r.TH2F(hname, ';p_{T} [GeV]; #eta', len(ptBinEdges)-1,  ptBinEdges, len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(v,
                  dict([(g,
                         dict([(s,
                                {'loose' : histo(var=v, hname=histoNamePerSamplePerSource(v, g, s, 'loose', mode)),
                                 'tight' : histo(var=v, hname=histoNamePerSamplePerSource(v, g, s, 'tight', mode)),
                                 })
                               for s in sources]))
                        for g in samples+['anygroup']]))
                 for v in variables])

def extractName(dictOrHist):
    "input must be either a dict or something with 'GetName'"
    isDict = type(dictOrHist) is dict
    return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()
def histoNamesPerSamplePerSource(variables, samples, leptonSources, mode) :
    return extractName(bookHistosPerSamplePerSource(variables, samples, leptonSources, mode=mode))

def writeHistos(outputFileName='', histosPerSamplePerSource={}, verbose=False):
    rootUtils.writeObjectsToFile(outputFileName, histosPerSamplePerSource, verbose)

def fetchHistos(fileName='', histoNames={}, verbose=False):
    return rootUtils.fetchObjectsFromFile(fileName, histoNames, verbose)

def computeEfficiency(histosTightLoose={}) :
    num = histosTightLoose['tight']
    den = histosTightLoose['loose']
    eff = num.Clone(num.GetName().replace('tight', 'tight_over_loose'))
    eff.Divide(den)
    return eff
def computeEfficiencies(dictOfHistosTightLoose={}):
    doh = dictOfHistosTightLoose
    isTlHistos = 'tight' in doh and 'loose' in doh
    return computeEfficiency(doh) if isTlHistos else dict((k, computeEfficiencies(v)) for k,v in doh.iteritems())

if __name__=='__main__':
    main()
