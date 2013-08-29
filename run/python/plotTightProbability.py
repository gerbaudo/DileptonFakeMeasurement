#!/bin/env python

# From the histograms produced by measureTightProbability, compute the
# conditional probability vs. pt and plot it.
#
# davide.gerbaudo@gmail.com
# Aug 2013

import collections
import glob
import optparse
import os
import re
import subprocess
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from rootUtils import referenceLine

#from datasets import datasets
from utils import getCommandOutput, guessMonthDayTag, longestCommonSubstring
from NavUtils import getAllHistoNames
import datasets

def buildRatioHisto(hNumerator, hDenominator) :
    hN, hD = hNumerator, hDenominator
    if not hD.GetEntries() : return
    # todo: rebin when low stats
    eff = r.TGraphAsymmErrors()
    lcs = longestCommonSubstring
    commonName  = lcs(hN.GetName(), hD.GetName())
    commonTitle = lcs(hN.GetTitle(), hD.GetTitle())
    eff.SetName('ptight_'+commonName)
    eff.SetTitle('p(tight) '+commonTitle)
    eff.Divide(hN, hD) # compute errors correctly, see TGraphAsymmErrors:Divide
    return eff
def findHistonamePairs(histonames, suffix1st='', suffix2nd='') :
    histonames = list(set(histonames))
    def name2nd(name1st, suff1st, suff2nd) : return name1st.replace(suff1st, suff2nd)
    pairs = [(f, s) for f, s in [(f, name2nd(f, suffix1st, suffix2nd))
                                 for f in histonames if f.endswith(suffix1st)]
             if s in histonames]
    return pairs
def buildProbHisto(file, hnameNum, hnameDen) :
    hNum, hDen = file.Get(hnameNum), file.Get(hnameDen)
    # hNum.Rebin(2)
    # hDen.Rebin(2)
    return buildRatioHisto(hNum, hDen)
def processFile(filename, outdir, label='') :
    file = r.TFile.Open(filename)
    outfname = (outdir+'/%(histoname)s_'
                +os.path.basename(filename).replace('.root','.png'))
    histonames = getAllHistoNames(file, onlyTH1=True)
    probHistos = [buildProbHisto(file, n, d)
                  for n,d in findHistonamePairs(histonames, '_num', '_den')]
    probHistos = filter(None, probHistos) # remove empty graphs that make TAxis complain
    c = r.TCanvas('p_tight','')
    for ph in probHistos :
        c.cd()
        c.Clear()
        if label : ph.SetTitle("%s %s"%(ph.GetTitle(), ", %s"%label if label else ''))
        ph.Draw('ap')
        xAx, yAx = ph.GetXaxis(), ph.GetYaxis()
        yAx.SetRangeUser(0.0, 1.1)
        l = referenceLine(xAx.GetXmin(), xAx.GetXmax())
        l.Draw()
        c.Update()
        c.SaveAs(outfname%{'histoname':ph.GetName()})

usage="""%prog [options] input1 [input2...]
Compute the conditional probability and plot it

Example:
.%prog -v out/fakeprob/merged/*<tag>.root
"""
parser = optparse.OptionParser(usage=usage)
parser.add_option('-o', '--outdir', help=('output directory;'
                                          +' default <input>/plotsTightProbability/'))
parser.add_option('-v', '--verbose', action='store_true', help='print details')
parser.add_option('-d', "--debug", action='store_true', help='print even more details')
(options, args) = parser.parse_args()
inputs, inputDir, ext = args, None, '.root'
if len(inputs) < 1 : parser.error("provide at least one input")
inputs = [f  for i in inputs for f in glob.glob(i if i.endswith(ext) else i+'/*'+ext)]
inputdir = os.path.dirname(inputs[0])
verbose  = options.verbose
debug    = options.debug
outdir   = options.outdir if options.outdir else inputdir+'/plotsTightProbability/'
tag      = guessMonthDayTag(inputs[0])

if verbose :
    print "Options:"
    print '\n'.join([str(inputs)]
                    + ["%s : %s" % (o, eval(o))
                       for o in ['inputdir', 'outdir', 'tag', 'verbose', 'debug',]])
if not os.path.isdir(outdir) :
    os.mkdir(outdir)
    if verbose : print "created directory '%s'"%outdir

allGroups = datasets.allGroups(datasets.datasets)
allSamples = datasets.allDatasets(datasets.datasets)
for f in inputs :
    group  = next((g for g in allGroups if g in f), None)
    sample = next((s for s in allSamples if s in f), None)
    label = "%s %s"%(group if group else '', sample if sample else '')
    processFile(f, outdir, label)

# allDatasets = [d for d in datasets if not d.placeholder]
# filenamesByGroup = collections.defaultdict(list)
# rootfiles = filter(os.path.isfile, glob.glob(inputdir + "*.root"))
# rootfiles = [rf for rf in rootfiles if tag in rf]

