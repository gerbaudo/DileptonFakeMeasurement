#!/bin/env python


import array
import glob
import os
import ROOT as r
r.gROOT.SetStyle('Plain')
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from utils import first
from rootUtils import drawLegendWithDictKeys
from SampleUtils import colors

rootcoredir = os.environ['ROOTCOREDIR']
r.gROOT.LoadMacro(rootcoredir+'/scripts/load_packages.C+')
r.load_packages()
tlv = r.TLorentzVector
#r.gInterpreter.GenerateDictionary('vector< vector< int > >', 'vector')

treename = 'SusySel'
basedir = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04_dev/SusyTest0/run/out/susysel/merged/'
samples = ['diboson', 'heavyflavor', 'ttbar', 'wjets', 'zjets']
filenames = dict((s, glob.glob(basedir+'/'+s+'*root')) for s in samples)
assert all(len(v)==1 for v in filenames.values()),"ambiguous filenames\n%s"%str(filenames)
filenames = dict((s,v[0]) for s,v in filenames.iteritems())

def FourMom2TLorentzVector(fm) :
    l = tlv()
    l.SetPxPyPzE(fm.px, fm.py, fm.pz, fm.E)
    return l

fm2tlv = FourMom2TLorentzVector
mt2 = r.mt2_bisect.mt2()
def computeMt2(a, b, met, zeroMass, lspMass) :
    pa    = array.array( 'd', [0.0 if zeroMass else a.M(), a.Px(), a.Py() ] )
    pb    = array.array( 'd', [0.0 if zeroMass else b.M(), b.Px(), b.Py() ] )
    pmiss = array.array( 'd', [0.0, met.Px(), met.Py() ] )
    mt2.set_momenta(pa, pb, pmiss)
    mt2.set_mn(lspMass)
    return mt2.get_mt2()
def computeMt2j(l0, l1, j0, j1, met, zeroMass=False, lspMass=0.0) :
    "As described in CMS-SUS-13-017"
    mt2_00 = computeMt2(l0+j0, l1+j1, met, zeroMass, lspMass)
    mt2_01 = computeMt2(l1+j0, l0+j1, met, zeroMass, lspMass)
    return min([mt2_00, mt2_01])
def computeMljj(l0, l1, j0, j1) :
    "Todo: extened combinatorics to N_j>2; good enough for now (we have few cases with >=3j)"
    jj = j0+j1
    dr0, dr1 = jj.DeltaR(l0), jj.DeltaR(l1)
    return (jj+l0).M() if dr0<dr1 else (jj+l1).M()

hs_mt2j = dict()
hs_mljj = dict()
for sample, filename in filenames.iteritems() :
    file = r.TFile.Open(filename)
    tree = file.Get(treename)
    print "processing %s (%d entries)"%(sample, tree.GetEntries())
    iEvent = 0
    h_mt2j = r.TH1F('h_mt2j_'+sample, ';m^{J}_{T2} [GeV]; entries/bin', 20, 0.0, 400.0)
    h_mljj = r.TH1F('h_mljj_'+sample, ';m_{ljj} [GeV]; entries/bin', 35, 0.0, 700.0)
    for event in tree :
        l0 = fm2tlv(event.l0)
        l1 = fm2tlv(event.l1)
        met = fm2tlv(event.met)
        jets = [fm2tlv(j) for j in event.jets]
        pars = event.pars
        if len(jets)<2 : continue
        j0, j1 = jets[0], jets[1]
        mt2_ja = computeMt2j(l0, l1, j0, j1, met)
        mt2_jb = computeMt2j(l0, l1, j0, j1, met, zeroMass=True)
        mt2_jc = computeMt2j(l0, l1, j0, j1, met, lspMass=100.0)
        mljj = computeMljj(l0, l1, j0, j1)
        weight, evtN = pars.weight, pars.eventNumber
        if iEvent<5 :
            print evtN,') mt2_j(m!=0): ',mt2_ja,', mt2_j(m==0): ',mt2_jb,', mt2_j(lsp=100): ',mt2_jc
        iEvent += 1
        h_mt2j.Fill(mt2_ja, weight)
        h_mljj.Fill(mljj, weight)
    hs_mt2j[sample] = h_mt2j
    hs_mljj[sample] = h_mljj


def plot(histos, var) :
    c = r.TCanvas('c_'+var,'')
    c.cd()
    pm = first(histos)
    pm.SetMaximum(1.1*max([h.GetMaximum() for h in histos.values()]))
    pm.Draw('axis')
    markers = dict(zip(histos.keys(), [r.kPlus, r.kCircle, r.kMultiply, r.kOpenSquare, r.kOpenTriangleUp]))
    for s,h in histos.iteritems() :
        h.SetLineColor(colors[s])
        h.SetMarkerColor(colors[s])
        h.SetLineWidth(2*h.GetLineWidth())
        h.SetMarkerStyle(markers[s])
        h.Draw('same')
    leg = drawLegendWithDictKeys(c, histos)
    leg.Draw()
    c.Update()
    c.SaveAs(var+'.png')

plot(hs_mt2j, 'mt2j')
plot(hs_mljj, 'mljj')
