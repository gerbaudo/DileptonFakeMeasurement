#!/bin/env python


# Compare the fake matrix histograms between two files produced by SusyTest0/FinalNewFake.
#
# davide.gerbaudo@gmail.com
# September 2013

import optparse
import os
import re
import ROOT as r

r.gROOT.SetBatch(1)
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
r.gErrorIgnoreLevel=r.kError               # disable root warnings 
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

from NavUtils import getAllHistoNames
from rootUtils import drawLegendWithDictKeys
from utils import filterWithRegexp

def main(filename1, filename2, outdir, regexp, verbose) :
    # there are too many histos; by default compare only the ones vs. pt for few selections
    relevantHistograms = ['l_pt_den','l_pt_num']
    relevantSelections = ['CRPremT2','CR_SSInc','CR_WHSS','CRPremT2','CR_SSInc','CR_WHSS']
    outdir = outdir if outdir else guessOutdirFromInputs(filename1, filename2)
    if verbose : print "saving output plots to '%s'"%outdir
    file1, file2 = r.TFile.Open(filename1), r.TFile.Open(filename2)
    histonames1 = getAllHistoNames(file1, onlyTH1=True)
    histonames2 = getAllHistoNames(file2, onlyTH1=True)
    if verbose :
        print '\n'.join(["%s: %d histograms"%(f, len(hs))
                         for f, hs in [(filename1, histonames1), (filename2, histonames2)]])
    commonHistos = [h for h in histonames1 if h in histonames2]
    def diff(a, b) : return list(set(a)-set(b))
    if len(commonHistos)!=len(histonames2) :
        missFrom1, missFrom2 = diff(histonames2, histonames1), diff(histonames1, histonames2)
        print ('histograms not in common:\n'
               +"%d missing from %s\n"%(len(missFrom1), filename1)
               #+'\n\t'.join(missFrom1)
               +"%d missing from %s\n"%(len(missFrom2), filename2)
               #+'\n\t'.join(missFrom2)
               )
    label1, label2 = labelFromFilename(filename1), labelFromFilename(filename2)
    commonHistos = (filterWithRegexp(commonHistos, regexp) if regexp is not None
                    else filter(lambda h:
                                any(s in h for s in relevantHistograms) and
                                any(s in h for s in relevantSelections),
                                commonHistos))
    canvas = r.TCanvas('diff_fakeMatrix','diff_fakeMatrix')
    for h in commonHistos :
        outname = outdir+'/'+h
        h1, h2 = file1.Get(h), file2.Get(h)
        plotComparison(h1, h2, canvas, outname, label1, label2, verbose)
def labelFromFilename(filename) :
    fname = os.path.basename(filename)
    fname, ext = os.path.splitext(fname)
    return fname
def plotComparison(histo1, histo2, canvas, outname, label1, label2, verbose) :
    h1, h2 = histo1, histo2
    if not h1.Integral() and not h2.Integral() : return
    canvas.cd()
    canvas.Clear()
    h1.SetMarkerColor(r.kBlack)
    h2.SetMarkerColor(r.kRed)
    h1.SetMarkerStyle(r.kFullCircle)
    h2.SetMarkerStyle(r.kOpenCircle)
    h1.SetMaximum(1.1*max([h.GetMaximum() for h in [h1, h2]]))
    h1.SetStats(0)
    h1.SetStats(0)
    h1.Draw()
    h2.Draw('same')
    ks = h1.KolmogorovTest(h2) if h2.Integral() and h1.Integral() else 0.0
    leg = drawLegendWithDictKeys(canvas, {label1 : h1, label2 : h2, },
                                 legWidth=0.35, legHeight=0.275)
    header = "%s : KS=%.3f"%(h1.GetName(), ks)
    leg.SetHeader(header)
    if ks<0.05 : print header
    canvas.Update()
    canvas.SaveAs(outname+'.png')
    if verbose : print "saving %s"%outname
def guessOutdirFromInputs(input1, input2) :
    def guessLocalInput(input1, input2) :
        return next(f for f in [input1, input2] if '../' not in f)
    return os.path.dirname(guessLocalInput(input1, input2))

if __name__=='__main__' :
    usage="""%prog input1 input2
    Plot a comparison and a ratio of the matrix histograms
    
    Example:
    %prog -v out/fakerate/merged/fakeout.root  ../../SusyMatrixMethod/data/forDavide_Sep11_2013.root
    """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-o', '--outdir', help=('output directory; default <input>/plots/'))
    parser.add_option('-r', '--regexp', default=None,
                      help=('only consider histo with names matching re, eg. \'pt\'.'
                            +'Default: only a subset of pt histos for few signal regions'))
    parser.add_option('-v', '--verbose', action='store_true', help='print details')
    (options, args) = parser.parse_args()
    inputs = args
    if len(inputs) < 2 : parser.error("provide two inputs")
    main(inputs[0], inputs[1], options.outdir, options.regexp, options.verbose)

