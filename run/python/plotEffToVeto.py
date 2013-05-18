#!/bin/env python

# plot, in the (mc1,mn1) plane the efficiency wrt. the base selection
#
# Inputs:
# - the pickle files with the number of signal events with and w/out veto
#   (produced with makeCutflowTable.py)
#
# davide.gerbaudo@gmail.com
# May 2013

import collections, optparse, re, sys
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from PickleUtils import readFromPickle
from SampleUtils import ModeAWhDbPar, ModeAWhDbReqid

#########
# default parameters [begin]
defaultSigVetoPickle   = 'counts_sig_May17_n0139_raw.pkl'
defaultSigNovePickle   = 'counts_sig_May17_n0139_raw_noveto.pkl'
defaultSigScale    = 1.0
# default parameters [end]
#########

parser = optparse.OptionParser()
parser.add_option("--sig-veto", dest="sigve", default=defaultSigVetoPickle,
                  help="file with signal counts with veto, default : %s" % defaultSigVetoPickle)
parser.add_option("--sig-nove", dest="signv", default=defaultSigNovePickle,
                  help="file with signal counts with veto, default : %s" % defaultSigNovePickle)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
sigInVetoFname  = options.sigve
sigInNoveFname  = options.signv
verbose         = options.verbose

countsSigVetoSampleSel = readFromPickle(sigInVetoFname)
countsSigNoveSampleSel = readFromPickle(sigInNoveFname)

reqDb = ModeAWhDbReqid()
parDb = ModeAWhDbPar()

def selIsRelevant(sel) : return any([sel.startswith(s) for s in ['sr6','sr7','sr8','sr9']])
def selIsBase(sel) : return 'base' in sel
def selIsFinal(sel) : return sel in ['sr6','sr7','sr8','sr9']
def getBaseSel(sel) : return re.search('(sr\d+)', sel).group(1)+'base'

def nicefySelectionName(s) :
    return s.replace('eq2j', ' N_{j}==2')\
           .replace('ge2j', ' N_{j}>=2')\
           .replace('ge3j', ' N_{j}>=3')\
           .replace('Nfv', ' no fw-jet veto')

mc1Range = {'min': min(parDb.allMc1()), 'max' : max(parDb.allMc1())}
mn1Range = {'min': min(parDb.allMn1()), 'max' : max(parDb.allMn1())}

histos = dict()
allNumeratorSelections = [s for s in list(set(k for countsSel in countsSigVetoSampleSel.values() for k in countsSel.keys()))
                          if selIsFinal(s)]
print allNumeratorSelections
allSignalSamples = countsSigVetoSampleSel.keys()

histos = dict()
for sel in allNumeratorSelections :
    histos[sel] = r.TH2F(sel+'_veto_nove',
                         nicefySelectionName(sel)+' : fraction yield veto(HTautau) / yield noveto. ;mc_{1};mn_{1}',
                         50, float(mc1Range['min']), float(mc1Range['max']),
                         50, float(mn1Range['min']), float(mn1Range['max']))

percent = 100.
#fill histo for signal points
for sample, countsSel in countsSigVetoSampleSel.iteritems() :
    mc1, mn1 = parDb.mc1Mn1ByReqid(reqDb.reqidBySample(sample))
    for sel, counts in sorted(countsSel.iteritems()) :
        if not selIsFinal(sel) : continue
        histo = histos[sel]
        refCounts = countsSigNoveSampleSel[sample][sel]
        if refCounts : histo.Fill(mc1, mn1, counts/refCounts)

# draw histos and print eff
r.gStyle.SetPaintTextFormat('.1f')
maxEff = 1.0
for s, h in histos.iteritems() :
    c = r.TCanvas('c_effVeto_'+s, 'relative eff '+s, 800, 600)
    c.cd()
    h.SetStats(0)
    h.SetMarkerSize(1.5*h.GetMarkerSize())
    h.SetMaximum(maxEff)
    h.Draw('colz') #contz
    h.SetMarkerSize(1.5*h.GetMarkerSize())
    h.Draw('text same')
    def writeSigMinMaxAvg(h, x=0.9, y=0.9) :
        tex = r.TLatex(0.0, 0.0, '')
        tex.SetNDC()
        tex.SetTextFont(h.GetTitleFont())
        tex.SetTextAlign(33)
        effs = [h.GetBinContent(i,j)
                for i in range(1,h.GetNbinsX()+1)
                for j in range(1,h.GetNbinsY()+1)
                if h.GetBinContent(i,j)]        
        tex.DrawLatex(x, y, 'sig : '+', '.join(["%s=%.1f"%(l,v) for l,v in [('min',min(effs)),
                                                                           ('avg',sum(effs)/len(effs) if len(effs) else 0.),
                                                                           ('max',max(effs))]]))
    writeSigMinMaxAvg(h)
    c.Update()
    c.SaveAs(c.GetName()+'.png')
