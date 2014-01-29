#!/bin/env python

# plot one single histogram from merged files produced with SusyPlot
#
# davide.gerbaudo@gmail.com
# Jan 2014

import datetime
from math import sqrt
import optparse
from rootUtils import (getNumDenHistos
                       ,getMinMax
                       ,topRightLegend
                       ,importRoot
                       )
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)
from utils import (first
                   ,rmIfExists
                   ,mkdirIfNeeded
                   )
from SampleUtils import colors

usage="""
Example usage:
%prog \\
 --histo histo_name \\
 --tag ${TAG} \\
 --input_dir  out/fakerate/merged/ \\
 --output_dir out/fakerate/merged/plot_${TAG} \\
 >& log/fakerate/FakePlot_${TAG}.log
"""
def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-n', '--histoname')
    parser.add_option('-o', '--output_dir')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_dir', 'histoname', 'output_dir']
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag            = opts.tag.strip('_')
    inputDirname   = opts.input_dir
    inputDirname   = inputDirname+'/' if not inputDirname.endswith('/') else inputDirname
    histoName      = opts.histoname
    outputDirname  = opts.output_dir
    outputDirname  = outputDirname+'/' if not outputDirname.endswith('/') else outputDirname
    mkdirIfNeeded(outputDirname)
    verbose        = opts.verbose
    if verbose : print ('\nUsing the following options:\n'
                        +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))

    inputFiles = getInputFiles(inputDirname, tag, verbose)
    assert all(f for f in inputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in inputFiles.iteritems()]))
    histograms = getHistograms(inputFiles, histoName)
    can = r.TCanvas('can_'+histoName, histoName, 800, 600)
    draw(can, histograms, label=histoName)
    can.SaveAs(outputDirname+'/'+histoName+'.png')

def mcSamples() : return ['ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def dataSample() : return 'data'
def susyplotSamples() : return [dataSample()] + mcSamples()
def getInputFiles(inputDirname, tag, verbose=False) :
    print "getInputFiles ~duplicated with buildWeightedMatrix.py; refactor"
    inDir = inputDirname
    samples = susyplotSamples()
    files = dict(zip(samples, [r.TFile.Open(inDir+'/'+s+'_'+tag+'.root') for s in samples]))
    if verbose : print "getInputFiles('%s'):\n\t%s"%(inputDirname, '\n\t'.join("%s : %s"%(k, f.GetName()) for k, f in files.iteritems()))
    return files
def getHistograms(inputFiles={}, histoName='') :
    "retrieve the nominal input histograms, build the tot hist, and compute the error"
    histograms = dict()
    for sample, file in inputFiles.iteritems() :
        h = file.Get(histoName)
        h.SetDirectory(0)
        assert h,"missing %s from %s"%(histoname, file.GetName())
        histograms[sample] = h
    return histograms

def draw(pad, hists, label=('','')) :
    pad.Draw()
    pad.cd()
    h_data = hists[dataSample()]
    stack = r.THStack('stack_'+h_data.GetName(), '')
    leg = topRightLegend(pad, 0.275, 0.475, shift=-0.025)
    leg.SetBorderSize(0)
    leg.AddEntry(h_data, dataSample(), 'P')
    for s in mcSamples() :
        h = hists[s]
        h.SetMarkerSize(0)
        h.SetFillColor(colors[s])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
        stack.Add(h)
    for s in mcSamples()[::-1] : leg.AddEntry(hists[s], s, 'F') # stack goes b-t, legend goes t-b
    h_data.SetMarkerStyle(r.kFullDotLarge)
    stack.Draw('hist')
    h_data.Draw('p same')
    leg.Draw('same')
    tex = r.TLatex()
    tex.SetNDC(True)
    label = "#splitline{%s}{%s}"%(label[0], label[1]) if len(label)==2 else label
    tex.DrawLatex(0.25, 0.95, label)
    pad.Update() # force stack to create padMaster
    padMaster = stack.GetHistogram()
    pMin, pMax = getMinMax([padMaster, h_data])
    stack.SetMaximum(pMax)
    pad._graphical_objects = [stack, h_data, leg, tex] + [h for h in stack.GetStack()]
    pad.Update()

if __name__=='__main__' :
    main()
