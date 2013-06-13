#!/bin/env python

# plot, in the (mc1,mn1) plane the samples' reqids
#
# Inputs:
# - the info stored in SampleUtils
#
# davide.gerbaudo@gmail.com
# June 2013

import collections, optparse, re, sys
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from PickleUtils import readFromPickle
from SampleUtils import ModeAWhDbPar, ModeAWhDbReqid, ModeAWhDbMergedFake2Lreqid


parser = optparse.OptionParser()
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
verbose         = options.verbose

reqDb = ModeAWhDbReqid()
parDb = ModeAWhDbPar()

allMc1 = [float(m) for m in parDb.allMc1()]
allMn1 = [float(m) for m in parDb.allMn1()]
def roundup(val) : return round(float(val)+0.5)
def rounddo(val) : return round(float(val)+0.0)
mc1Range = {'min': rounddo(min(allMc1)), 'max' : roundup(max(allMc1))}
mn1Range = {'min': rounddo(min(allMn1)), 'max' : roundup(max(allMn1))}

histo2l = r.TH2F('mc1mn1_2lep_reqids',
                 'ReqIds for the WH 2lep grid ;mc_{1};mn_{1}',
                 50, float(mc1Range['min']), float(mc1Range['max']),
                 50, float(mn1Range['min']), float(mn1Range['max']))
histo3l = r.TH2F('mc1mn1_3lep_reqids',
                 'ReqIds for the WH 3lep grid ;mc_{1};mn_{1}',
                 50, float(mc1Range['min']), float(mc1Range['max']),
                 50, float(mn1Range['min']), float(mn1Range['max']))
# # looping over reqDb.entries also includes missing entries (i.e. reqids that were not processed)
#  for sample, reqid in reqDb.entries.iteritems() :
#      mc1, mn1 = parDb.mc1Mn1ByReqid(reqid)
for entry in parDb.entries :
    reqid, mc1, mn1 = entry.ds, entry.mc1, entry.mn1
    mc1, mn1 = float(mc1), float(mn1)
    sample = reqDb.sampleByReqid(reqid)
    if not sample :
        print "missing %s, (%.1f, %.1f)"%(reqid, mc1, mn1)
        continue
    h = histo2l if '2Lep' in sample else histo3l if '3Lep' in sample else None
    assert h,"%s is not 2l nor 3l"%sample
    h.Fill(mc1, mn1, float(reqid))

def buildMergedHisto(histo) :
    merge = ModeAWhDbMergedFake2Lreqid()
    h = histo.Clone(histo.GetName()+'_merged')
    h.Reset()
    xAx, yAx = h.GetXaxis(), h.GetYaxis()
    for iX in range(1, h.GetNbinsX()+1) :
        for iY in range(1, h.GetNbinsY()+1) :
            x,y = xAx.GetBinCenter(iX), yAx.GetBinCenter(iY)
            fakeReqid = merge.reqidByMc1Mn1(x, y)
            if fakeReqid : h.Fill(x, y, fakeReqid)
    return h

r.gStyle.SetPaintTextFormat('.0f')
c = r.TCanvas('c_reqids', 'WH reqids', 800, 600)
c.cd()
for h in [histo2l, histo3l] :
    c.Clear()
    h.SetTitle(h.GetTitle()+" (%d points)"%h.GetEntries())
    h.SetStats(0)
    h.SetMarkerSize(1.25*h.GetMarkerSize())
    hMerge = buildMergedHisto(h)
    hMerge.Draw('col')
    h.Draw('text20 same')
    c.Update()
    c.SaveAs('c_'+h.GetName()+'.png')
