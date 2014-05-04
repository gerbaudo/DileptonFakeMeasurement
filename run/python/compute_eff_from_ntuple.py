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
    optionsToPrint = ['inputDir', 'outputDir', 'mode', 'tag', 'doFillHistograms', 'cacheFileName']
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
    groups = [g for g in samplesPerGroup.keys() if not isDataSample(g)]
    sourcesThisMode = {'real' : ['real'], # use same convention as in FakeLeptonSources.h
                       'conv' : ['conv'],
                       'hflf' : ['heavy', 'light', 'qcd']
                       }[mode]
    #fill histos
    if doFillHistograms :
        histosPerGroupPerSource = bookHistosPerSamplePerSource(vars, groups, sourcesThisMode, mode=mode)
        for group in groups:
            filenames = filenamesPerGroup[group]
            histosThisGroupPerSource = dict((v, histosPerGroupPerSource[v][group]) for v in histosPerGroupPerSource.keys())
            chain = r.TChain(treeName)
            [chain.Add(fn) for fn in filenames]
            if verbose: print "%s : %d entries"%(group, chain.GetEntries())
            fillHistos(chain, histosThisGroupPerSource, lepton, mode, verbose)
        writeHistos(cacheFileName, histosPerGroupPerSource, verbose)
    # compute efficiencies
    histosPerGroupPerSource = fetchHistos(cacheFileName, histoNamesPerSamplePerSource(vars, groups, sourcesThisMode, mode), verbose)
    effs = computeEfficiencies(histosPerGroupPerSource) # still [var][gr][source][l/t]
    for s in sourcesThisMode:
        for v in vars:
            varIs1D = v=='pt'
            if varIs1D:
                effsThisSourceThisVar = dict((g, effs[v][g][s]) for g in groups)
                cname = 'eff_'+lepton+'_'+s
                lT, lX, lY = '#varepsilon(T|L)', 'p_{T} [GeV]', '#varepsilon(T|L)'
                title = lT+' '+s+' '+lepton+';'+lX+';'+lY
                zoomIn = True
                fakeu.plot1dEfficiencies(effsThisSourceThisVar, cname, outputDir, title, zoomIn)
    writeHistos(outputFileName, effs, verbose)
    if verbose : print "saved scale factors to %s" % outputFileName

#___________________________________________________

leptonTypes = ['loose', 'tight']
leptonSources = []
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def fillHistos(chain, histosPerSource, lepton, mode, verbose=False):
    class Counters: # scope trick (otherwise unavailable within nested func
        nLoose, nTight = 0, 0
        totWeightLoose, totWeightTight = 0.0, 0.0
        def str(self):
            counterNames = ['nLoose', 'nTight', 'totWeightLoose', 'totWeightTight']
            return ', '.join(["%s : %.1f"%(c, getattr(self, c)) for c in counterNames])
    counters = Counters()
    addTlv = kin.addTlv
    for iEvent, event in enumerate(chain) :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber       
        l0, l1 = event.l0, event.l1
        def fillHistosBySource(lep):
            isTight = lep.isTight
            source = enum2source(lep)
            isRightLep = lep.isEl if lepton=='el' else probe.isMu
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
                if mode=='hflf': # qcd is heavy+light
                    histosPerSource['pt'    ]['qcd'][tightOrLoose].Fill(pt,      weight)
                    histosPerSource['pt_eta']['qcd'][tightOrLoose].Fill(pt, eta, weight)
            if isRightLep and isRightSource :
                lep = addTlv(lep)
                pt, eta = lep.p4.Pt(), abs(lep.p4.Eta())
                fill('loose')
                if isTight : fill('tight')
        fillHistosBySource(l0)
        fillHistosBySource(l1)
    if verbose : print counters.str()

def histoNamePerSamplePerSource(var, sample, leptonSource, tightOrLoose, mode):
    return 'h_'+var+'_'+sample+'_'+leptonSource+'_'+tightOrLoose+'_'+mode
def bookHistosPerSamplePerSource(variables, samples, sources, mode=''):
    "book a dict of histograms with keys [var][sample][lepton_source][tight, loose]"
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
                        for g in samples]))
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

def plotStackedHistos(histosPerGroup={}, outputDir='', mode='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    leptonTypes = first(first(histosPerGroup)).keys()
    colors = SampleUtils.colors
    histosPerName = dict([(mode+'_'+var+'_'+lt, # one canvas for each histo, so key with histoname w/out group
                           dict([(g, histosPerGroup[g][var][lt]) for g in groups]))
                          for var in variables for lt in leptonTypes])
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
        topRightLabel(can, histoname, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]+[data]
        can.Update()
        can.SaveAs(os.path.join(outputDir, histoname+'.png'))
def plotStackedHistosSources(histosPerVar={}, outputDir='', mode='', verbose=False):
    variables = histosPerVar.keys()
    sources = first(histosPerVar).keys()
    colors = colorsFillSources
    for var in variables:
        histos = dict((s, histosPerVar[var][s]['loose']) for s in sources)
        canvasBasename = mode+'_region_'+var+'_loose'
        missingSources = [s for s, h in histos.iteritems() if not h]
        if missingSources:
            if verbose : print "skip %s, missing histos for %s"%(var, str(missingSources))
            continue
        totBkg = summedHisto(histos.values())
        err_band = buildErrBandGraph(totBkg, computeStatErr2(totBkg))
        emptyBkg = totBkg.Integral()==0
        if emptyBkg:
            if verbose : print "empty backgrounds, skip %s"%canvasBasename
            continue
        can = r.TCanvas('c_'+canvasBasename, canvasBasename, 800, 600)
        can.cd()
        pm = totBkg # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update() # necessary to fool root's dumb object ownership
        stack = r.THStack('stack_'+canvasBasename,'')
        can.Update()
        r.SetOwnership(stack, False)
        for s, h in histos.iteritems() :
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
        topRightLabel(can, canvasBasename, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, dictSum(histos, {'stat err':err_band}), opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]
        can.Update()
        can.SaveAs(os.path.join(outputDir, canvasBasename+'.png'))

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
def plotPerSourceEff(histosPerVar={}, outputDir='', lepton='', mode='', sample='', verbose=False, zoomIn=True):
    "plot efficiency for each source (and 'anysource') as a function of each var; expect histos[var][source][loose,tight]"
    variables = histosPerVar.keys()
    sources = [s for s in first(histosPerVar).keys() if s!='real'] # only fake eff really need a scale factor
    colors = colorsLineSources
    for var in filter(lambda x : x in ['pt1', 'eta1'], histosPerVar.keys()):
        histosPerSource = dict((s, histosPerVar[var][s]) for s in sources)
        canvasBasename = mode+'_efficiency_'+lepton+'_'+var+("_%s"%sample if sample else '')
        can = r.TCanvas('c_'+canvasBasename, canvasBasename, 800, 600)
        can.cd()
        pm = first(effs) # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update()
        for s, h in effs.iteritems() :
            h.SetMarkerColor(colors[s] if s in colors else r.kBlack)
            h.SetLineColor(h.GetMarkerColor())
            h.SetLineWidth(2*h.GetLineWidth())
            h.SetMarkerStyle(markersSources[s] if s in markersSources else r.kDot)
            h.Draw('ep same')
            h.SetDirectory(0)
        #pprint.pprint(effs)
        yMin, yMax = getMinMax(effs.values())
        pm.SetMinimum(0.0)
        pm.SetMaximum(0.25 if yMax < 0.5 and zoomIn else 1.1)
        can.Update()
        topRightLabel(can, canvasBasename, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, effs, opt='lp')
        can.RedrawAxis()
        can._histos = effs
        can.Update()
        can.SaveAs(os.path.join(outputDir, canvasBasename+'.png'))

def subtractRealAndComputeScaleFactor(histosPerGroup={}, variable='', outhistoname='', outputDir='./', mode='', verbose=False):
    "efficiency scale factor"
    groups = histosPerGroup.keys()
    histosPerType = dict([(lt,
                           dict([(g,
                                  histosPerGroup[g][variable][lt])
                                 for g in groups]))
                          for lt in leptonTypes])
    for lt in leptonTypes :
        histosPerType[lt]['totSimBkg'] = summedHisto([histo for group,histo in histosPerType[lt].iteritems() if isBkgSample(group)])

    simuTight = histosPerType['fake_tight']['totSimBkg']
    simuLoose = histosPerType['fake_loose']['totSimBkg']
    dataTight = histosPerType['tight'     ]['data'     ]
    dataLoose = histosPerType['loose'     ]['data'     ]
    # subtract real contribution from data
    # _Note to self_: currently estimating the real contr from MC; in
    # the past also used iterative corr, which might be more
    # appropriate in cases like here, where the normalization is
    # so-so.  Todo: investigate the normalization.
    dataTight.Add(histosPerType['real_tight']['totSimBkg'], -1.0)
    dataLoose.Add(histosPerType['real_loose']['totSimBkg'], -1.0)
    dataTight.Divide(dataLoose)
    simuTight.Divide(simuLoose)
    print "eff(T|L) vs. ",variable
    def formatFloat(floats): return ["%.4f"%f for f in floats]
    print "efficiency data : ",formatFloat(getBinContents(dataTight))
    print "efficiency simu : ",formatFloat(getBinContents(simuTight))
    ratio = dataTight.Clone(outhistoname)
    ratio.SetDirectory(0)
    ratio.Divide(simuTight)
    print "sf    data/simu : ",formatFloat(getBinContents(ratio))
    print "            +/- : ",formatFloat(getBinErrors(ratio))
    can = r.TCanvas('c_'+outhistoname, outhistoname, 800, 600)
    botPad, topPad = rootUtils.buildBotTopPads(can)
    can.cd()
    topPad.Draw()
    topPad.cd()
    pm = dataTight
    pm.SetStats(0)
    pm.Draw('axis')
    xAx, yAx = pm.GetXaxis(), pm.GetYaxis()
    xAx.SetTitle('')
    xAx.SetLabelSize(0)
    yAx.SetRangeUser(0.0, 0.25)
    textScaleUp = 1.0/topPad.GetHNDC()
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetTitle('#epsilon(T|L)')
    yAx.SetTitleOffset(yAx.GetTitleOffset()/textScaleUp)
    simuTight.SetLineColor(r.kRed)
    simuTight.SetMarkerStyle(r.kOpenCross)
    simuTight.SetMarkerColor(simuTight.GetLineColor())
    dataTight.Draw('same')
    simuTight.Draw('same')
    leg = drawLegendWithDictKeys(topPad, {'data':dataTight, 'simulation':simuTight}, legWidth=0.4)
    leg.SetHeader('scale factor '+mode+' '+('electron' if '_el_'in outhistoname else 'muon' if '_mu_' in outhistoname else ''))
    can.cd()
    botPad.Draw()
    botPad.cd()
    ratio.SetStats(0)
    ratio.Draw()
    textScaleUp = 1.0/botPad.GetHNDC()
    xAx, yAx = ratio.GetXaxis(), ratio.GetYaxis()
    yAx.SetRangeUser(0.0, 2.0)
    xAx.SetTitle({'pt1':'p_{T}', 'eta1':'|#eta|'}[variable])
    yAx.SetNdivisions(-202)
    yAx.SetTitle('Data/Sim')
    yAx.CenterTitle()
    xAx.SetLabelSize(textScaleUp*0.04)
    xAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    refLine = rootUtils.referenceLine(xAx.GetXmin(), xAx.GetXmax())
    refLine.Draw()
    can.Update()
    can.SaveAs(os.path.join(outputDir, mode+'_'+outhistoname+'.png'))
    can.SaveAs(os.path.join(outputDir, mode+'_'+outhistoname+'.eps'))
    return ratio

def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
#     print "computeStatErr2 use the one in rootUtils"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}

def buildErrBandGraph(histo_tot_bkg, err2s) :
#     print "buildErrBandGraph use the one in rootUtils"
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
