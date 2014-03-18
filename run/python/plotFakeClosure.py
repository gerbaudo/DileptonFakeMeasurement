#!/bin/env python

# Plot histograms for fake closure test

# Inputs: histograms from root files from SusyPlot and from FakePred

# add description here
# ref to Matt's code

# davide.gerbaudo@gmail.com
# October 2013

import numpy as np
import optparse
import os
from rootUtils import (referenceLine
                       ,topRightLegend
                       ,getMinMax
                       ,buildRatioHistogram
                       ,buildBotTopPads
                       ,importRoot
                       ,integralAndError
                       )
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
from utils import (enumFromHeader
                   ,json_write
                   ,mkdirIfNeeded
                   ,rmIfExists
                   )
from SampleUtils import colors
from systUtils import (buildErrBandGraph
                       ,buildErrBandRatioGraph
                       ,computeStatErr2
                       ,fetchFakeSysHistosAndComputeSysErr2
                       )


usage="""
Example usage:
%prog \\
 --tag ${TAG} \\
 --input_dir  out/susyplot/merged/ \\
 --input_fake out/fakepred/merged/data_${TAG}.root \\
 --output_dir out/fakepred/merged/fake_closure_plots_${TAG} \\
 >& log/fakerate/FinalFakeHist_${TAG}.log
"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-f', '--input_fake')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-o', '--output_dir')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_fake', 'input_dir', 'output_dir',]
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag           = opts.tag.strip('_')
    inputFakeFile = opts.input_fake
    inputDirname  = opts.input_dir
    outputDir     = opts.output_dir
    outputDir     = outputDir if outputDir.endswith('/') else outputDir+'/'
    verbose       = opts.verbose
    if verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    inputFiles = getInputFiles(inputDirname, tag, verbose)
    inputFiles[fakeSample()] = r.TFile.Open(inputFakeFile)
    assert all(f for f in inputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in inputFiles.iteritems()]))
    mkdirIfNeeded(outputDir)

    for region in [#'cr8lptee', 'cr8lptmm', 'cr9lpt', 'crSsEwkLoose'
                   #,"crZVfake1jee"
                   #,"crZVfake2jee"
                   #,"crZVfake1jem"
                   #,"crZVfake2jem"
                   #,"crfake1jem"
                   #,"crfake2jem"
                   #,"crZV1jmm"
                   #,"crZV2jmm"
                   #,"crfake1jmm"
                   #,"crfake2jmm"

                   "crZVfake1j"
                   ,"crZVfake2j"
                   ,"crfake1j"
                   ,"crfake2j"
                   ,"crZV1j"
                   ,"crZV2j"

                   ] :
        for channel in ['ee', 'em', 'mm'] :
            for varname in ['l0_pt', 'l1_pt', 'll_M', 'metrel', 'met', 'njets', 'nbjets'] :
                histo_basename = region+'_'+channel+'_'+varname
                hists, err2s = buildHists(inputFiles, histo_basename)
                if not hists[dataSample()].GetEntries() : continue
                err_band     = buildErrBandGraph(hists['sm'], err2s)
                err_band_r   = buildErrBandRatioGraph(err_band)
                can = r.TCanvas('can_'+histo_basename, histo_basename, 800, 600)
                botPad, topPad = buildBotTopPads(can)
                can.cd()
                topPad.Draw()
                drawTop(topPad, hists, err_band, (channel, region))
                can.cd()
                botPad.Draw()
                drawBot(botPad, hists[dataSample()], hists['sm'], err_band_r, xaxisLabel(varname))
                can.Update()
                outFilename = outputDir+histo_basename+'.png'
                rmIfExists(outFilename) # avoid root warnings
                can.SaveAs(outFilename)
    if verbose : print "output saved to \n%s"%outputDir

def mcSamples() : return ['ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def bkSamples() : return ['fake']+mcSamples()
def dataSample() : return 'data'
def fakeSample() : return 'fake'
def isFake(samplename) : return fakeSample()==samplename
def isMc(samplename) : return samplename in mcSamples()
def isBkgSample(samplename) : return isFake(samplename) or isMc(samplename)
def susyplotSamples() : return [dataSample()] + mcSamples()
def getInputFiles(inputDirname, tag, verbose=False) :
    print "getInputFiles ~duplicated with buildWeightedMatrix.py; refactor"
    inDir = inputDirname
    samples = susyplotSamples()
    files = dict(zip(samples, [r.TFile.Open(inDir+'/'+s+'_'+tag+'.root') for s in samples]))
    if verbose : print "getInputFiles('%s'):\n\t%s"%(inputDirname, '\n\t'.join("%s : %s"%(k, f.GetName()) for k, f in files.iteritems()))
    return files
def xaxisLabel(varname) :
    return {'l0_pt'   : 'l_{0} P_{T} [GeV]'
            ,'l1_pt'  : 'l_{1} P_{T} [GeV]'
            ,'ll_M'   : 'm(ll) [GeV]'
            ,'met'    : '#slash{E}_{T} [GeV]'
            ,'metrel' : '#slash{E}^{rel}_{T} [GeV]'
            ,'njets'  : '# jets'
            ,'nbjets' : '# b jets'
            }[varname]
def cloneAndReset(h, clone_name='') :
    c = h.Clone(clone_name)
    c.SetDirectory(0)
    c.Reset()
    return c
def buildHists(inputFiles={}, histo_basename='') :
    "retrieve the nominal input histograms, build the tot hist, and compute the error"
    histograms = dict()
    h_sm, h_data = None, None
    err2s = {}
    for sample, file in inputFiles.iteritems() :
        histoname = histo_basename+('_NONE' if isFake(sample) else '_NOM')
        h = file.Get(histoname)
        assert h,"cannot get %s from %s"%(histoname, file.GetName())
        h.SetDirectory(0)
        if not h : print "=> missing %s from %s"%(histoname, file.GetName())
        assert h,"missing %s from %s"%(histoname, file.GetName())
        histograms[sample] = h
        if isFake(sample) : err2s = addErr2s(err2s, fetchFakeSysHistosAndComputeSysErr2(file, h))
        elif isMc(sample) : err2s = addErr2s(err2s, computeStatErr2(h))
        h_sm = h_sm if h_sm else cloneAndReset(h, 'SM') # just pick the first one as a template
        if isBkgSample(sample) : h_sm.Add(h, 1.0)
    histograms['sm'] = h_sm
    return histograms, err2s

def addErr2s(current_err2={}, add_err2={}) :
    if not current_err2 : return add_err2
    else :
        return dict([(k, current_err2[k] + add_err2[k]) for k in current_err2.keys()])

def drawTop(pad, hists, err_band, label=('','')) :
    pad.Draw()
    pad.cd()
    h_data = hists[dataSample()]
    h_bkg  = hists['sm']
    stack = r.THStack('stack_'+h_data.GetName(), '')
    leg = topRightLegend(pad, 0.275, 0.475, shift=-0.025)
    leg.SetBorderSize(0)
    leg.AddEntry(h_data, dataSample(), 'P')
    leg.AddEntry(h_bkg, 'sm', 'L')
    counters = dict()
    if 'nbjets' in h_data.GetName() :
        counters['data'] = integralAndError(h_data)
        counters['sm'] = integralAndError(h_bkg)
    for s in bkSamples() :
        h = hists[s]
        h.SetMarkerSize(0)
        h.SetFillColor(colors[s])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
        if 'nbjets' in h.GetName() : counters[s] = integralAndError(h)
        stack.Add(h)
    for s in bkSamples()[::-1] : leg.AddEntry(hists[s], s, 'F') # stack goes b-t, legend goes t-b
    h_data.SetMarkerStyle(r.kFullDotLarge)
    stack.Draw('hist')
    h_bkg.Draw('hist same')
    err_band.Draw('E2 same')
    h_data.Draw('p same')
    leg.AddEntry(err_band, 'Uncertainty', 'F')
    leg.Draw('same')
    tex = r.TLatex()
    tex.SetNDC(True)
    label = "#splitline{%s}{%s}"%(label[0], label[1]) if len(label)==2 else label
    tex.DrawLatex(0.50, 0.85, label)
    pad.Update() # force stack to create padMaster
    padMaster = stack.GetHistogram()
    pMin, pMax = getMinMax([h_bkg, h_data, err_band])
    stack.SetMaximum(pMax)
    xAx, yAx = padMaster.GetXaxis(), padMaster.GetYaxis()
    xAx.SetTitle('')
    xAx.SetLabelSize(0)
    textScaleUp = 1.0/pad.GetHNDC()
    yAx.SetLabelSize(textScaleUp*yAx.GetLabelSize())
    pad._graphical_objects = [stack, h_data, h_bkg, err_band, leg, tex] + [h for h in stack.GetStack()]
    pad.Update()
    if 'nbjets' in h_data.GetName() :
        print h_data.GetName()
        print " | ".join(["%12s"%k for k,v in sorted(counters.items())])
        print " | ".join(["%12s"%("%.1f +/- %.1f"%v) for k,v in sorted(counters.items())])


def getXrange(h) :
    nbins = h.GetNbinsX()
    x_lo = h.GetBinCenter(1) - 0.5*h.GetBinWidth(1)
    x_hi = h.GetBinCenter(nbins) + 0.5*h.GetBinWidth(nbins)
    return x_lo, x_hi
def drawBot(pad, histo_data, histo_mc, err_band_r, xaxis_title='') :
    pad.Draw()
    pad.cd()
    ratio = buildRatioHistogram(histo_data, histo_mc)
    yMin, yMax = 0.0, 2.0
    ratio.SetMinimum(yMin)
    ratio.SetMaximum(yMax)
    ratio.SetStats(0)
    ratio.Draw('axis')
    x_lo, x_hi = getXrange(ratio)
    refLines = [referenceLine(x_lo, x_hi, y, y) for y in [0.5, 1.0, 1.5]]
    for l in refLines : l.Draw()
    err_band_r.Draw('E2 same')
    ratio.Draw('ep same')
    xA, yA = ratio.GetXaxis(), ratio.GetYaxis()
    textScaleUp = 1.0/pad.GetHNDC()
    if xaxis_title : xA.SetTitle(xaxis_title)
    yA.SetNdivisions(-104)
    yA.SetTitle('Data/SM')
    yA.CenterTitle()
    yA.SetTitleOffset(yA.GetTitleOffset()/textScaleUp)
    for a in [xA, yA] :
        a.SetLabelSize(a.GetLabelSize()*textScaleUp)
        a.SetTitleSize(a.GetTitleSize()*textScaleUp)
    pad._graphical_objects = [ratio, err_band_r] + refLines # avoid garbage collection
    pad.Update()


if __name__=='__main__' :
    main()
