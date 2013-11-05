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
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
from rootUtils import (referenceLine
                       ,topRightLegend
                       ,getMinMax
                       )
from utils import (enumFromHeader
                   ,json_write
                   ,rmIfExists
                   )
from SampleUtils import colors

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
    tag = opts.tag
    inputFakeFile = opts.input_fake
    inputDirname  = opts.input_dir
    outputDir     = opts.output_dir
    verbose       = opts.verbose
    if verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    inputFiles = getInputFiles(inputDirname, tag, verbose)
    inputFiles[fakeSample()] = r.TFile.Open(inputFakeFile)
    assert all(f for f in inputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in inputFiles.iteritems()]))

    for region in ['cr8lptee', 'cr8lptmm', 'cr9lpt', 'cr8lptmmMtww', 'cr8lptmmHt'] :
        for channel in ['ee', 'em', 'mm'] :
            for varname in ['l0_pt', 'l1_pt', 'll_M', 'metrel', 'met', 'njets', 'nbjets'] :
                xlabel = xaxisLabel(varname)
                histo_basename = region+'_'+channel+'_'+varname
                hists, err2s = buildHists(inputFiles, histo_basename)
                err_band     = buildErrBandGraph(hists['sm'], err2s)
                err_band_r   = buildErrBandRatioGraph(err_band)
                can = r.TCanvas('can_'+histo_basename, histo_basename, 800, 600)
                botPad, topPad = buildBotTopPads(can)
                can.cd()
                topPad.Draw()
                drawTop(topPad, hists, err_band, (channel, region))
                can.cd()
                botPad.Draw()
                drawBot(botPad, hists[dataSample()], hists['sm'], err_band_r)
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
    tag = tag if tag.startswith('_') else '_'+tag
    samples = susyplotSamples()
    files = dict(zip(samples, [r.TFile.Open(inDir+'/'+s+tag+'.root') for s in samples]))
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
        h.SetDirectory(0)
        if not h : print "=> missing %s from %s"%(histoname, file.GetName())
        assert h,"missing %s from %s"%(histoname, file.GetName())
        histograms[sample] = h
        if isFake(sample) : err2s = addErr2s(err2s, computeFakeSysErr2(file, h))
        elif isMc(sample) : err2s = addErr2s(err2s, computeStatErr2(h))
        h_sm = h_sm if h_sm else cloneAndReset(h, 'SM') # just pick the first one as a template
        if isBkgSample(sample) : h_sm.Add(h, 1.0)
    histograms['sm'] = h_sm
    return histograms, err2s

def addErr2s(current_err2={}, add_err2={}) :
    if not current_err2 : return add_err2
    else :
        return dict([(k, current_err2[k] + add_err2[k]) for k in current_err2.keys()])

def computeFakeSysErr2(input_fake_file=None, nominal_histo=None) :
    "Compute the bin-by-bin sum2 err including up&down fake systematic variations + stat. unc."
    variations = ['_'+l+'_'+s+'_'+ud # 2x2x2=8 syst variations for the fake estimate
                  for l in ['EL', 'MU'] for s in ['FR', 'RE'] for ud in ['UP','DOWN']]
    variations = ['_EL_RE_UP', '_EL_RE_DOWN', '_MU_RE_UP', '_MU_RE_DOWN', '_EL_FR_UP',
                  '_EL_FR_DOWN', '_MU_FR_UP', '_MU_FR_DOWN']
    nom_hname = nominal_histo.GetName()
    vars_histos = dict([(v, input_fake_file.Get(nom_hname.replace('_NONE',v))) for v in variations])
    def bc(h) : return [h.GetBinContent(b) for b in range(1, 1+h.GetNbinsX())]
    def be(h) : return [h.GetBinError(b) for b in range(1, 1+h.GetNbinsX())]
    nom_bcs  = bc(nominal_histo)
    nom_be2s = [e*e for e in be(nominal_histo)]
    vars_bcs = dict([(v, bc(h)) for v, h in vars_histos.iteritems()])
    bins = range(nominal_histo.GetNbinsX())
    deltas = [[vars_bcs[v][b] - nom_bcs[b] for v in vars_bcs.keys()] for b in bins]
    def positive(ll) : return [l if not l<0.0 else 0.0 for l in ll ]
    def negative(ll) : return [l if     l<0.0 else 0.0 for l in ll ]
    def sumquad(ll) : return sum([l*l for l in ll])
    up_e2s = np.array([sumquad(positive(deltas[b])) for b in bins]) + np.array(nom_be2s)
    do_e2s = np.array([sumquad(negative(deltas[b])) for b in bins]) + np.array(nom_be2s)
    return {'up' : up_e2s, 'down' : do_e2s}
def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}

def buildErrBandGraph(histo_tot_bkg, err2s) :
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
def buildErrBandRatioGraph(errband_graph) :
    gr = errband_graph.Clone()
    points = range(gr.GetN())
    xs     = np.array([gr.GetX()[i] for i in points])
    ys     = np.array([gr.GetY()[i] for i in points])
    eys_lo = np.array([abs(gr.GetErrorYlow (i)) for i in points])
    eys_hi = np.array([abs(gr.GetErrorYhigh(i)) for i in points])
    ys_lo  = ys - eys_lo
    ys_hi  = ys + eys_hi
    def absFracDeltaFromUnity(y_nom, y_var) : return abs(y_var/y_nom - 1.0) if y_nom else 0.0
    eys_lo = [absFracDeltaFromUnity(n, v) for n, v in zip(ys, ys_lo)]
    eys_hi = [absFracDeltaFromUnity(n, v) for n, v in zip(ys, ys_hi)]
    for p, x, ey_lo, ey_hi in zip(points, xs, eys_lo, eys_hi) :
        gr.SetPoint(p, x, 1.0) # TGraph does not have a SetPointY, so we need to set both x and y
        gr.SetPointEYlow (p, ey_lo)
        gr.SetPointEYhigh(p, ey_hi)
    return gr
def buildBotTopPads(canvas, splitFraction=0.275) :
    canvas.cd()
    botPad = r.TPad(canvas.GetName()+'_bot', 'bot pad', 0.0, 0.0, 1.0, splitFraction, 0, 0, 0)
    interPadMargin = 0.5*0.05
    botPad.SetTopMargin(interPadMargin)
    botPad.SetBottomMargin(botPad.GetBottomMargin()/splitFraction)
    botPad.SetRightMargin(0.20*botPad.GetRightMargin())
    r.SetOwnership(botPad, False)
    canvas.cd()
    canvas.Update()
    topPad = r.TPad(canvas.GetName()+'_top', 'top pad', 0.0, splitFraction, 1.0, 1.0, 0, 0)
    topPad.SetBottomMargin(interPadMargin)
    topPad.SetTopMargin(0.20*topPad.GetTopMargin())
    topPad.SetRightMargin(0.20*topPad.GetRightMargin())
    r.SetOwnership(topPad, False)
    canvas._pads = [topPad, botPad]
    return botPad, topPad

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
    for s in bkSamples() :
        h = hists[s]
        h.SetMarkerSize(0)
        h.SetFillColor(colors[s])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
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
    tex.DrawLatex(0.45, 0.75, label)
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

def buildRatioHistogram(num, den) :
    print "buildRatioHistogram: move to utils"
    ratio = num.Clone(num.GetName()+'_over_'+den.GetName())
    ratio.SetDirectory(0) # we usually don't care about the ownership of these temporary objects
    ratio.Reset()
    ratio.Divide(num, den, 1, 1, 'B')
    return ratio
def getXrange(h) :
    nbins = h.GetNbinsX()
    x_lo = h.GetBinCenter(1) - 0.5*h.GetBinWidth(1)
    x_hi = h.GetBinCenter(nbins) + 0.5*h.GetBinWidth(nbins)
    return x_lo, x_hi
def drawBot(pad, histo_data, histo_mc, err_band_r) :
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
