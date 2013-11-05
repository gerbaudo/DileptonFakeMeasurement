#!/bin/env python

# Root utility functions
#
# davide.gerbaudo@gmail.com
# 2013-08-26


import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
import numpy as np
from utils import verticalSlice

def referenceLine(xmin=0., xmax=100.0, ymin=1.0, ymax=1.0) :
    l1 = r.TLine(xmin, ymin, xmax, ymax)
    l1.SetLineStyle(3)
    l1.SetLineColor(r.kGray+1)
    return l1
def firstHisto(histos) :
    return (histos.itervalues().next() if type(histos) is dict
            else histos[0] if type(histos) is list
            else None)
def unitLineFromFirstHisto(histos) :
    fH = firstHisto(histos)
    xAx = fH.GetXaxis()
    return referenceLine(xAx.GetXmin(), xAx.GetXmax())

def topRightLegend(pad,  legWidth, legHeight, shift=0.0) :
    rMarg, lMarg, tMarg = pad.GetRightMargin(), pad.GetLeftMargin(), pad.GetTopMargin()
    leg = r.TLegend(1.0 - rMarg - legWidth + shift,
                    1.0 - tMarg - legHeight + shift,
                    1.0 - rMarg + shift,
                    1.0 - tMarg + shift)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    pad._leg = leg
    return leg
def drawLegendWithDictKeys(pad, histosDict, legWidth=0.325, legHeight=0.225) :
    leg = topRightLegend(pad, legWidth, legHeight)
    for s,h in histosDict.iteritems() :
        leg.AddEntry(h, s, 'p')
    leg.Draw()
    pad.Update()
    return leg
def getMinMaxFromTGraph(gr) :
    points = range(gr.GetN())
    y = np.array([gr.GetY()[p] for p in points])
    return min(y), max(y)
def getMinMaxFromTGraphAsymmErrors(gr) :
    points = range(gr.GetN())
    y    = np.array([gr.GetY()[p] for p in points])
    y_el = np.array([abs(gr.GetErrorYlow (i)) for i in points])
    y_eh = np.array([abs(gr.GetErrorYhigh(i)) for i in points])
    return min(y-y_el), max(y+y_eh)
def getMinMaxFromTH1(h) :
    bins = range(1, 1+h.GetNbinsX())
    y   = np.array([h.GetBinContent(b) for b in bins])
    y_e = np.array([h.GetBinError(b)   for b in bins])
    return min(y-y_e), max(y+y_e)
def getMinMax(histosOrGraphs=[]) :
    def mM(obj) :
        cname = obj.Class().GetName()
        if   cname.startswith('TH1') :               return getMinMaxFromTH1(obj)
        elif cname.startswith('TGraphAsymmErrors') : return getMinMaxFromTGraphAsymmErrors(obj)
        elif cname.startswith('TGraph') :            return getMinMaxFromTGraph(obj)
    ms, Ms = verticalSlice([mM(o) for o in histosOrGraphs])
    return min(ms), max(Ms)
