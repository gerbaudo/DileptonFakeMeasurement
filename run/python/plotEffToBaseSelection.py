#!/bin/env python

# plot, in the (mc1,mn1) plane the efficiency wrt. the base selection
#
# Inputs:
# - the pickle files with the number of signal and background events
#   (produced with makeCutflowTable.py)
#
# davide.gerbaudo@gmail.com
# Feb 2013

import collections, optparse, re, sys
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from PickleUtils import readFromPickle
from SampleUtils import ModeAWhDbPar, ModeAWhDbReqid

#########
# default parameters [begin]
defaultSigPickle   = 'counts_signal.pkl'
defaultBkgPickle   = 'counts_backgr.pkl'
defaultSigScale    = 1.0
# default parameters [end]
#########

parser = optparse.OptionParser()
parser.add_option("-s", "--sig-file", dest="sig", default=defaultSigPickle,
                  help="file with signal counts, default : %s" % defaultSigPickle)
parser.add_option("-b", "--bkg-file", dest="bkg", default=defaultBkgPickle,
                  help="file with background counts, default : %s" % defaultBkgPickle)
parser.add_option("-c", "--counts", action='store_true', dest="cnt", default=False,
                  help="plot counts rather than efficiencies")
parser.add_option("-r", "--raw-counts", action='store_true', dest="rawcnt", default=False,
                  help="plot raw counts rather than efficiencies")
parser.add_option("-S", "--scale-sig", dest="sigScale", default=defaultSigScale, type='float',
                  help="scale the signal yield by this factor (default %.1f)" % defaultSigScale)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
sigInputFname   = options.sig
bkgInputFname   = options.bkg
cnt             = options.cnt
rawcnt          = options.rawcnt
sigScale        = options.sigScale
verbose         = options.verbose

doCounts = cnt or rawcnt
doEff = not doCounts

countsSigSampleSel = readFromPickle(sigInputFname)
countsBkgSampleSel = readFromPickle(bkgInputFname)

countBkgTot = collections.defaultdict(float)
for sample, countsSel in countsBkgSampleSel.iteritems() :
    for sel, counts in countsSel.iteritems() :
        countBkgTot[sel] += counts
countsBkgSampleSel['totBkg'] = countBkgTot

reqDb = ModeAWhDbReqid()
parDb = ModeAWhDbPar()

def selIsRelevant(sel) : return any([sel.startswith(s) for s in ['sr6','sr7','sr8','sr9']])
def selIsBase(sel) : return 'base' in sel
def getBaseSel(sel) : return re.search('(sr\d+)', sel).group(1)+'base'

def nicefySelectionName(s) :
    return s.replace('eq2j', ' N_{j}==2')\
           .replace('ge2j', ' N_{j}>=2')\
           .replace('ge3j', ' N_{j}>=3')\
           .replace('Nfv', ' no fw-jet veto')

mc1Range = {'min': min(parDb.allMc1()), 'max' : max(parDb.allMc1())}
mn1Range = {'min': min(parDb.allMn1()), 'max' : max(parDb.allMn1())}

histos = dict()
allNumeratorSelections = [s for s in list(set(k for k in countsSel.keys() for countsSel in countsSigSampleSel.values()))
                          if selIsRelevant(s) and not selIsBase(s)]
print allNumeratorSelections
allSignalSamples = countsSigSampleSel.keys()

histos = dict()

for sel in allNumeratorSelections :
    if not selIsRelevant(sel) or selIsBase(sel) : continue
    baseSel = getBaseSel(sel)
    title = nicefySelectionName(sel)
    title += (' : eff. [%] rel. to '+baseSel+';mc_{1};mn_{1}') if  doEff else ''
    title += ' signal counts;mc_{1};mn_{1}' if doCounts and cnt else ''
    title += ' signal raw counts;mc_{1};mn_{1}' if doCounts and rawcnt else ''
    histos[sel] = r.TH2F(sel+'_'+baseSel, title,
                         50, float(mc1Range['min']), float(mc1Range['max']),
                         50, float(mn1Range['min']), float(mn1Range['max']))

percent = 100.
#fill histo for signal points
for sample, countsSel in countsSigSampleSel.iteritems() :
    mc1, mn1 = parDb.mc1Mn1ByReqid(reqDb.reqidBySample(sample))
    print "%s (%.1f, %.1f) " % (sample, mc1, mn1)
    for sel, counts in sorted(countsSel.iteritems()) :
        if not selIsRelevant(sel) or selIsBase(sel) : continue
        baseSel = getBaseSel(sel)
        histo = histos[sel]
        refCounts = countsSigSampleSel[sample][baseSel]
        if doCounts            : histo.Fill(mc1, mn1, counts)
        if doEff and refCounts : histo.Fill(mc1, mn1, percent*counts/refCounts)

# draw histos and print bkg eff
r.gStyle.SetPaintTextFormat('.0f' if rawcnt else '.2f')
maxEff = 100. if doEff else max([h.GetMaximum() for h in histos.values()])
for s, h in histos.iteritems() :
    c = r.TCanvas(('c_eff_' if doEff else 'c_counts_raw_' if rawcnt else 'c_counts_')+s, ''+s, 800, 600)
    c.cd()
    h.SetStats(0)
    h.SetMarkerSize(1.5*h.GetMarkerSize())
    h.SetMaximum(maxEff)
    h.Draw('colz') #contz
    h.SetMarkerSize(1.5*h.GetMarkerSize())
    h.Draw('text same')
    def writeBkgEff(backgrounds, counts,
                    font='', x=0.125, y=0.9, slope=-0.05) :
        tex = r.TLatex(0.0, 0.0, '')
        tex.SetNDC()
        if font : tex.SetTextFont(font)
        tex.SetTextAlign(13)
        for bkg, cnts in zip(backgrounds,counts) :
            
            #tex.DrawLatex(1.0-c.GetTopMargin(), 1.0-c.GetRightMargin(), "Signal: x %.1f"%scale)
            selCnt, refCnt =  cnts[s], cnts[getBaseSel(s)]
            tex.DrawLatex(x, y, bkg +': '+ ('-' if not refCnt else "%.1f%%"%(percent*selCnt/refCnt)))
            y += slope
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
    backgrounds = countsBkgSampleSel.keys()
    counts = [countsBkgSampleSel[b] for b in backgrounds]
    if doEff : writeBkgEff(backgrounds, counts, font=h.GetTitleFont())
    #writeSigMinMaxAvg(h)
    c.Update()
    c.SaveAs(c.GetName()+'.png')

