#!/usr/bin/env python

# Given an input file.root, loop over the histograms and rename them
# swapping two substrings (eg. old-label new-label); save to out.root
#
# davide.gerbaudo@gmail.com
# Oct 2013

import os
import sys
import ROOT as r
r.gROOT.SetBatch(True)
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from NavUtils import getAllHistoNames

def main(filenameOrig, filenameDest, substr1, substr2) :
    s1, s2 = substr1, substr2
    input = r.TFile.Open(filenameOrig)
    histonames = getAllHistoNames(input)
    output = r.TFile.Open(filenameDest, 'recreate')
    output.cd()
    nhFixed = 0
    for hn in histonames :
        h = input.Get(hn)
        hnn = hn.replace(s1,s2) if s1 in hn else hn.replace(s2,s1) # avoid swapping twice s1->s2->s1
        h.Write(hnn)
        if hn!=hnn : nhFixed += 1
    print "renamed %d histograms"%nhFixed
    output.Close()
    input.Close()


if __name__=='__main__' :
    if len(sys.argv)<5 :
        print """Usage:
        \t %s in.root out.root string1 string2
        Example:
        \t %sout/fakepred/merged/data_Oct_08.root out/fakepred/merged/fixed/data_Oct_08.root sr8_ sr8lpt_ls
        """ % sys.argv[0]
    filenameOrig = sys.argv[1]
    filenameDest = sys.argv[2]
    substr1      = sys.argv[3]
    substr2      = sys.argv[4]
    main(filenameOrig, filenameDest, substr1, substr2)
