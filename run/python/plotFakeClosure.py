#!/bin/env python

# Plot histograms for fake closure test

# Inputs: 

# davide.gerbaudo@gmail.com
# October 2013

import optparse
import os
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from utils import (enumFromHeader
                   ,json_write
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
    tag = opts.tag
    inputFakeFile = opts.input_fake
    inputDirname  = opts.input_dir
    outputDir     = opts.output_dir
    verbose       = opts.verbose
    if verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    inputFiles = getInputFiles(inputDirname, tag, verbose)
    inputFiles[fakeSample()] = r.TFile.Open(inputFakeFile)
    assert all(f for f in allInputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in allInputFiles.iteritems()]))

    for region in ['cr8lptmm', 'cr8lptee', 'cr9lpt', 'cr8lptmmMtww', 'cr8lptmmHt'] :
        for channel in ['ee', 'em', 'mm'] :
            for plotname in ['l0_pt', 'l1_pt', 'll_M', 'metrel', 'met', 'njets', 'nbjets'] :
                xlabel = xaxisLabel(plotname)
                hists = buildHists()
                # top pad
                # bot pad

    if verbose : print "output saved to \n%s"%outputDir

def mcSamples() : return ['ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def bkSamples() : return ['fake']+mcSamples()
def dataSample() : return 'data'
def fakeSample() : return 'fake'
def susyplotSamples() : return [dataSample()] + mcSamples()
def getInputFiles(inputDirname, tag, verbose=False) :
    print "getInputFiles ~duplicated with buildWeightedMatrix.py; refactor"
    inDir = inputDirname
    tag = tag if tag.startswith('_') else '_'+tag
    samples = susyplotSamples()
    files = dict(zip(samples, [r.TFile.Open(inDir+'/'+s+tag+'.root') for s in samples]))
    if verbose : print "getInputFiles('%s'):\n\t%s"%(inputDirname, '\n\t'.join("%s : %s"%(k, f.GetName()) for k, f in files.iteritems()))
    return files
def xaxisLabel(plotname) :
    return {'l0_pt'   : 'l_{0} P_{T} [GeV]'
            ,'l1_pt'  : 'l_{1} P_{T} [GeV]'
            ,'ll_M'   : 'm(ll) [GeV]'
            ,'met'    : '#slash{E}_{T} [GeV]'
            ,'metrel' : '#slash{E}^{rel}_{T} [GeV]'
            ,'njets'  : '# jets'
            ,'nbjets' : '# b jets'
            }[plotname]

if __name__=='__main__' :
    main()
