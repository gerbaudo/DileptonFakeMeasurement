#!/bin/env python

# script to compute the fake factor for the emu selection
# the fake factor is described in section C of ATL-COM-PHYS-2012-1166

# davide.gerbaudo@gmail.com
# Sept 2014

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
from compute_fake_el_scale_factor import buildErrBandGraph, computeStatErr2

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
    parser.add_option('-m', '--mode', help='emu')
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
    validModes = ['emu']
    if mode not in validModes : parser.error("invalid mode %s"%mode)
    tupleStem, treeName = filter(lambda _: _[0]==mode, fakeu.tupleStemsAndNames)[0]

    templateInputFilename = "*_%(stem)s_tuple_%(tag)s.root" % {'tag':tag, 'stem':tupleStem}
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
        print 'input filenames: ',os.path.join(inputDir, templateInputFilename)
    # collect inputs
    tupleFilenames = glob.glob(os.path.join(inputDir, templateInputFilename))
    samples = setSameGroupForAllData(fastSamplesFromFilenames(tupleFilenames, verbose))
    samplesPerGroup = collections.defaultdict(list)
    filenamesPerGroup = collections.defaultdict(list)
    mkdirIfNeeded(outputDir)
    for s, f in zip(samples, tupleFilenames) :
        samplesPerGroup[s.group].append(s)
        filenamesPerGroup[s.group].append(f)
    vars = ['pt', 'pt_eta']
    groups = [g for g in samplesPerGroup.keys() if g is not 'higgs']
    if lepton=='el' : groups = [g for g in groups if g is not 'heavyflavor']
    sourcesThisMode = ['real', 'conv', 'heavy', 'light', 'unknown'] if lepton=='el' else ['real', 'heavy', 'light', 'unknown']
    #fill histos
    if doFillHistograms :
        start_time = time.clock()
        num_processed_entries = 0
        histosPerGroupPerSource = bookHistosPerSamplePerSource(vars, groups, sourcesThisMode, mode=mode)
        for group in groups:
            filenames = filenamesPerGroup[group]
            sources = histosPerGroupPerSource.keys()
            histosThisGroupPerSource = dict((s, histosPerGroupPerSource[s][group]) for s in sources)
            histosAnyGroupPerSource  = dict((s, histosPerGroupPerSource[s]['anygroup']) for s in sources) if group!='data' else {}

            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            if verbose: print "%s : %d entries"%(group, chain.GetEntries())
            is_data = group in ['data']
            print 'is_data ',is_data
            num_processed_entries += fillHistos(chain=chain,
                                                histosPerSource=histosThisGroupPerSource,
                                                histosPerSourceAnygroup=histosAnyGroupPerSource,
                                                lepton=lepton,
                                                onthefly_tight_def=onthefly_tight_def,
                                                verbose=verbose)
        writeHistos(cacheFileName, histosPerGroupPerSource, verbose)
        end_time = time.clock()
        delta_time = end_time - start_time
        one_minute = 60
        if verbose:
            print ("processed {0:d} entries ".format(num_processed_entries)
                   +"in "+("{0:d} min ".format(int(delta_time/60)) if delta_time>one_minute else
                           "{0:.1f} s ".format(delta_time))
                   +"({0:.1f} kHz)".format(num_processed_entries/delta_time))
    # plot histos
    histosPerGroupPerSource = fetchHistos(cacheFileName, histoNamesPerSamplePerSource(vars, groups, sourcesThisMode, mode), verbose)

    # effs = computeEfficiencies(histosPerGroupPerSource) # still [var][gr][source][l/t]
    for v in vars:
        varIs1D, varIs2D = v=='pt', v=='pt_eta'
        densThisSourceThisVar = dictSum(dict((s, histosPerGroupPerSource[v]['anygroup'][s]['loose']) for s in sourcesThisMode),
                                        {'data' : histosPerGroupPerSource[v]['data']['unknown']['loose']})
        numsThisSourceThisVar = dictSum(dict((s, histosPerGroupPerSource[v]['anygroup'][s]['tight']) for s in sourcesThisMode),
                                        {'data' : histosPerGroupPerSource[v]['data']['unknown']['tight']})
        if varIs1D:
            lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#varepsilon(T|L)'
            cname = 'stack_loose_'+lepton
            lT, lY = 'loose '+lepton+', denominator to #varepsilon(T|L)', '#varepsilon(T|L)'
            title = lT+' '+'anysource'+' '+lepton+';'+lX+';'+lY
            plotStackedHistosWithData(densThisSourceThisVar,
                                      outputDir, cname, title,
                                      colors=fakeu.colorsFillSources(),
                                      verbose=verbose)
            cname = 'stack_tight_'+lepton
            lT, lY = 'tight '+lepton+', numerator to #varepsilon(T|L)', '#varepsilon(T|L)'
            title = lT+' '+'anysource'+' '+lepton+';'+lX+';'+lY
            plotStackedHistosWithData(numsThisSourceThisVar,
                                      outputDir, cname, title,
                                      colors=fakeu.colorsFillSources(),
                                      verbose=verbose)

    for s in sourcesThisMode:
        for v in vars:
            groups = first(histosPerGroupPerSource).keys()
            varIs1D, varIs2D = v=='pt', v=='pt_eta'
            # effsThisSourceThisVar = dict((g, effs[v][g][s]) for g in groups)
            densThisSourceThisVar = dictSum(dict((g, histosPerGroupPerSource[v][g][s]['loose'])
                                                 for g in groups if g not in ['anygroup','data']),
                                            {'data' : histosPerGroupPerSource[v]['data']['unknown']['loose']})
            numsThisSourceThisVar = dictSum(dict((g, histosPerGroupPerSource[v][g]['unknown']['tight'])
                                                 for g in groups if g not in ['anygroup','data']),
                                            {'data' : histosPerGroupPerSource[v]['data']['unknown']['tight']})
            if varIs1D:
                # cname = 'eff_'+lepton+'_'+s
                lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#varepsilon(T|L)'
                # title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                # zoomIn = True
                # fakeu.plot1dEfficiencies(effsThisSourceThisVar, cname, outputDir, title, zoomIn)
                cname = 'stack_loose_'+lepton+'_'+s
                lT, lY = 'loose '+lepton+', denominator to #varepsilon(T|L)', '#varepsilon(T|L)'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                plotStackedHistosWithData(densThisSourceThisVar,
                                          outputDir, cname, title,
                                          colors=SampleUtils.colors,
                                          verbose=verbose)
                cname = 'stack_tight_'+lepton+'_'+s
                lT, lY = 'tight '+lepton+', numerator to #varepsilon(T|L)', '#varepsilon(T|L)'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                plotStackedHistosWithData(numsThisSourceThisVar,
                                          outputDir, cname, title,
                                          colors=SampleUtils.colors,
                                          verbose=verbose)

            # elif varIs2D:
            #     cname = 'eff_'+lepton+'_'+s
            #     lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#eta'
            #     title = lT+' '+s+' '+lepton+';'+lX+';'+lY
            #     fakeu.plot2dEfficiencies(effsThisSourceThisVar, cname, outputDir, title, zoomIn=zoomIn)
    # writeHistos(outputFileName, effs, verbose)
    if verbose : print "saved scale factors to %s" % outputFileName

#___________________________________________________

leptonTypes = ['loose', 'tight']
leptonSources = []
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def fillHistos(chain, histosPerSource, histosPerSourceAnygroup={}, lepton='', onthefly_tight_def=None, verbose=False):
    """fill the histograms, returns the number of events
    processed. histosPerSource is required; histosPerSourceAnygroup is
    filled only when provided"""
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
    sources = histosPerSource['pt'].keys() # tmp debug
    for iEvent, event in enumerate(chain) :
        num_processed_entries += 1
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        l0, l1 = event.l0, event.l1
        def is_tight(lep) : return onthefly_tight_def(lep) if onthefly_tight_def else lep.isTight
        def fillHistosBySource(lep):
            lep = addTlv(lep)
            isTight = is_tight(lep)
            source = enum2source(lep)
            isRightLep = lep.isEl if lepton=='el' else lep.isMu
            def fill(tightOrLoose):
                if tightOrLoose=='loose':
                    counters.nLoose +=1
                    counters.totWeightLoose += weight
                if tightOrLoose=='tight':
                    counters.nTight +=1
                    counters.totWeightTight += weight
                histosPerSource['pt'    ][source][tightOrLoose].Fill(pt,      weight)
                histosPerSource['pt_eta'][source][tightOrLoose].Fill(pt, eta, weight)
                if hspsAnygroup:
                    hspsAnygroup['pt'    ][source][tightOrLoose].Fill(pt,      weight)
                    hspsAnygroup['pt_eta'][source][tightOrLoose].Fill(pt, eta, weight)
            if isRightLep:
                lep = addTlv(lep)
                pt, eta = lep.p4.Pt(), abs(lep.p4.Eta())
                fill('loose')
                if isTight : fill('tight')
        isEmu = (l0.isEl and l1.isMu) or (l0.isMu and l1.isEl)
        l0, l1 = addTlv(l0), addTlv(l1)
        l0, l1 = (l0, l1) if l0.p4.Pt()>l1.p4.Pt() else (l1, l0) # pt sort
        # todo: add trigger requirement on l0
        isTightLoose = is_tight(l0) and not is_tight(l1)
        isTightTight = (is_tight(l0) and is_tight(l1))
        if isEmu and (isTightTight or isTightLoose) : fillHistosBySource(l1)
    if verbose : print counters.str()
    return num_processed_entries

def histoNamePerSamplePerSource(var, sample, leptonSource, tightOrLoose, mode):
    return 'h_'+var+'_'+sample+'_'+leptonSource+'_'+tightOrLoose+'_'+mode
def bookHistosPerSamplePerSource(variables, samples, sources, mode='', verbose=False):
    "book a dict of histograms with keys [var][sample][lepton_source][tight, loose]"
    if verbose : print 'booking histos'
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

def plotStackedHistosWithData(histosPerGroup={}, outputDir='', canvasname='', canvastitle='', colors={}, verbose=False):
    "histosPerGroup[group], where group=data is treated as special"
    groups = histosPerGroup.keys()
    mkdirIfNeeded(outputDir)
    missingGroups = [g for g, h in histosPerGroup.iteritems() if not h]
    if missingGroups:
        if verbose : print "skip %s, missing histos for %s"%(histoname, str(missingGroups))
        return
    bkgHistos = dict([(g, h) for g, h in histosPerGroup.iteritems() if not isDataSample(g)])
    totBkg = summedHisto(bkgHistos.values())
    err_band = buildErrBandGraph(totBkg, computeStatErr2(totBkg))
    emptyBkg = totBkg.Integral()==0
    histoname, region = totBkg.GetName(), 'emu' # tmp replacement vars, to be fixed
    if emptyBkg:
        if verbose : print "empty backgrounds, skip %s"%histoname
        return
    can = r.TCanvas(canvasname, canvastitle, 800, 600)
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
        if verbose :
            print "integrals : {0} tot.bkg.: {1}, data: {2}".format(histoname, totBkg.Integral(), data.Integral())
    else:
        print "no data"
    yMin, yMax = getMinMax([h for h in [totBkg, data, err_band] if h])
    pm.SetMinimum(0.0)
    pm.SetMaximum(1.1*yMax)
    can.Update()
    topRightLabel(can, "#splitline{%s}{%s}"%(histoname, region), xpos=0.15, ypos=(1.0-0.5*can.GetTopMargin()), align=13)
    drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
    can.RedrawAxis()
    can._stack = stack
    can._histos = [h for h in stack.GetHists()]+[data]
    can.Update()
    filename=os.path.join(outputDir, histoname+'.png')
    rmIfExists(filename)
    can.SaveAs(filename)




if __name__=='__main__':
    main()
