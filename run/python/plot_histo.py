#!/bin/env python


# given a histo name, and input files, plot the stack
# if only filenames are given, print the list of available histos
# the sample name is guessed from the filename (should be something coming out of mergeOutput.py)
#
# davide.gerbaudo@gmail.com
# Mar 2014

import collections
import glob
import optparse
import os
from pprint import pformat
import re
import subprocess
from datasets import datasets, setSameGroupForAllData
from rootUtils import (referenceLine
                       ,topRightLegend
                       ,getMinMax
                       ,buildRatioHistogram
                       ,buildBotTopPads
                       ,importRoot
                       )
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
from utils import (getCommandOutput,
                   first,
                   )
from NavUtils import getAllHistoNames
import SampleUtils

stackedBackgrounds = ['diboson', 'higgs', 'wjets', 'ttbar', 'zjets', 'heavyflavor']
backgrounds = stackedBackgrounds + ['allBkgButHf', 'allBkg']
samples = backgrounds + ['data']
colors = dict(SampleUtils.colors.items() + {'allBkgButHf' : r.kGray, 'allBkg' : r.kGray, 'higgs' : r.kMagenta}.items())

def main() :

    parser = optparse.OptionParser()
    parser.add_option("-n", "--histo-name")
    parser.add_option("-o", "--output-dir", default='./')
    parser.add_option("-s", "--suffix")
    (options, args) = parser.parse_args()
    histoName = options.histo_name
    outputDir = options.output_dir
    suffix    = options.suffix
    inputFileNames = args

    print '\n'.join(inputFileNames)
    sampleForFilename = dict((f, guessSampleFromFilename(f, samples)) for f in inputFileNames)
    assert all(s for f,s in sampleForFilename.iteritems())," cannot identify sample for some inputs %s"%pformat(sampleForFilename)
    filenames = collections.defaultdict(list)
    for fn, s in sampleForFilename.iteritems() : filenames[s].append(fn)
    assert all(len(fnames)==1 for s,fnames in filenames.iteritems()),"expected one file per sample, got %s"%pformat(filenames)
    filenames = dict((s, first(fn)) for s, fn in filenames.iteritems())
    files = dict((s, r.TFile.Open(fn)) for s, fn in filenames.iteritems())

    if histoName :
        print "plotting %s"%histoName
        canvasname = histoName+("_%s"%suffix if suffix else '')
        plotHisto(files, histoName, canvasname, outputDir)
    else :
        printHistoNames(first(files))

def guessSampleFromFilename(filename='', samples=[]) :
    "Guess sample from filename, slightly different from the one in SampleUtils.py"
    def longest(lst=[]) : return max(lst, key=len) if lst else None
    def filenameMatchesSample(fn, s) : return s in fn
    return longest([s for s in samples if filenameMatchesSample(filename, s)])

def plotHisto(inputFiles={}, histoname='', canvasname='', outdir='./') :
    can = r.TCanvas(canvasname, histoname, 800, 600)
    can.cd()
    backds = [b for b in backgrounds if b in inputFiles]
    stackBkgds = filter(lambda x : x in backds, stackedBackgrounds)
    h_bkgs = dict((b, inputFiles[b].Get(histoname)) for b in backds)
    h_btot = inputFiles['allBkg'].Get(histoname) if 'allBkg' in inputFiles else None
    h_bnoHf = inputFiles['allBkgButHf'].Get(histoname) if 'allBkgButHf' in inputFiles else None
    h_data = inputFiles['data'].Get(histoname) if 'data' in inputFiles else None
    stack = r.THStack('stack_'+histoname, '')
    leg = topRightLegend(can, 0.275, 0.475, shift=-0.025)
    if h_data : leg.AddEntry(h_data, "data: %.1f"%h_data.Integral(), 'P')
    if h_btot : leg.AddEntry(h_btot, "tot bkg %.1f"%h_btot.Integral(), 'l')
    if h_bnoHf : leg.AddEntry(h_bnoHf, "bkgNoHf: %.1f"%h_bnoHf.Integral(), 'l')
    for b in stackBkgds :
        h = h_bkgs[b]
        if not h :
            print 'missing h for ',b
            continue
        h.SetMarkerSize(0)
        h.SetFillColor(colors[b])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
        stack.Add(h)
    for s in stackBkgds[::-1] : leg.AddEntry(h_bkgs[s], "%s: %.1f"%(s, h_bkgs[s].Integral()), 'F') # stack goes b-t, legend goes t-b
    stack.Draw('hist')
    if h_btot  :
        h_btot.SetLineWidth(2*h_btot.GetLineWidth())
        h_btot.Draw('hist same')
    if h_bnoHf :
        h_bnoHf.SetLineWidth(2*h_bnoHf.GetLineWidth())
        h_bnoHf.SetLineStyle(2)
        h_bnoHf.Draw('hist same')
    if h_data  :
        h_data.SetMarkerStyle(r.kFullDotLarge)
        h_data.Draw('p same')
    leg.Draw('same')
    yMin, yMax = getMinMax([h for h in [stack.GetHistogram(), h_btot, h_bnoHf, h_data] if h is not None])
    stack.SetMaximum(1.1*yMax)
    tex = r.TLatex()
    tex.SetNDC(True)
    tex.DrawLatex(0.15, 0.95, histoname)
    can.Update()
    for ext in ['png'] : can.SaveAs(outdir+'/'+canvasname+'.'+ext)


def printHistoNames(inputFile) :
    print 'histograms from ',inputFile.GetName()
    print '\n'.join(getAllHistoNames(inputFile))
    
if __name__=='__main__' :
    main()
