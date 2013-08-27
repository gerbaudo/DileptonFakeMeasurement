#!/bin/env python

# Root utility functions
#
# davide.gerbaudo@gmail.com
# 2013-08-26


import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options

def referenceLine(xmin=0., xmax=100.0) :
    y1, y2 = 1.0, 1.0
    l1 = r.TLine(xmin, y1, xmax, y2)
    l1.SetLineStyle(3)
    l1.SetLineColor(r.kGray+1)
    return l1
