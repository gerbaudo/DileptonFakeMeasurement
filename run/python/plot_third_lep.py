#!/bin/env python

# plot the third lepton that we want to use to reject the diboson bkg
#
# Input:
# small ntuples from susy::wh::TupleMaker
# Output:
# plots of multiplicity and best candidate for a Z->ll lepton
#
# todo: include functions from utils, rootUtils, and avoid duplication
#
# davide.gerbaudo@gmail.com
# 2013-12-09

import array
import glob
import math
import os
import ROOT as r
r.gROOT.SetStyle('Plain')
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from kin import phi_mpi_pi, addTlv, lepIsSeparatedFromOther, lepPairIsZcand, deltaMZ0

r.gROOT.LoadMacro('src/TupleMakerObjects.h+')

treename = 'SusySel'
tag = 'Dec_09'
basedir = 'data/'

colors = {
    'ttbar'       : r.kRed+1,
    'zjets'       : r.kOrange-2,
    'wjets'       : r.kBlue-2,
    'diboson'     : r.kSpring+2,
    'singletop'   : r.kAzure-4,
    'multijet'    : r.kGray,
    'fake'        : r.kGray, # just another name for the same thing
    'heavyflavor' : r.kViolet+1
    }

def getInputFiles() :
    samples = ['diboson', 'WH_2Lep_3']
    filenames = dict((s, glob.glob(basedir+'/*'+s+'*'+tag+'.root')) for s in samples)
    assert all(len(v)==1 for v in filenames.values()),"ambiguous filenames\n%s"%str(filenames)
    filenames = dict((s,v[0]) for s,v in filenames.iteritems())
    return filenames

def first(listOrDict) :
    lod = listOrDict
    return lod.itervalues().next() if type(lod) is dict else lod[0] if lod else None


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
def drawLegendWithDictKeys(pad, histosDict, legWidth=0.325, legHeight=0.225, opt='p') :
    leg = topRightLegend(pad, legWidth, legHeight)
    for s,h in histosDict.iteritems() :
        leg.AddEntry(h, s, opt)
    leg.Draw()
    pad.Update()
    return leg
def fillHistos() :
    filenames = getInputFiles()
    hs_nlep = dict()
    hs_nillep = dict()
    hs_npairs = dict()
    hs_masses = dict()
    hs_bestm = dict()
    hs_bestlpt = dict()
    hs_bestleta= dict()
    hs_bestdphi= dict()
    hs_bestdrj = dict()
    for sample, filename in filenames.iteritems() :
        file = r.TFile.Open(filename)
        tree = file.Get(treename)
        print "processing %s (%d entries)"%(sample, tree.GetEntries())
        iEvent = 0
        h_nlep   = r.TH1F('h_nlep_'+sample, '; N_{lep, low-p_{T}};entries/bin', 5, -0.5, 4.5)
        h_nillep = r.TH1F('h_nillep_'+sample,
                          'separated from signal leptons; N_{lep, iso-lep, low-p_{T}};entries/bin',
                          5, -0.5, 4.5)
        h_npairs  = r.TH1F('h_npairs_'+sample, '; N_{valid pairs}; entries/bin', 5, -0.5, 4.5)
        h_masses  = r.TH1F('h_masses_'+sample, '; m_{ll}, any valid pair [GeV]; entries/bin', 25, 0.0, 200.0)
        h_bestm   = r.TH1F('h_bestm_'+sample, '; m_{ll}, best pair [GeV]; entries/bin', 25, 0.0, 200.0)
        h_bestlpt = r.TH1F('h_bestlpt_'+sample, '; p_{T}, best soft lepton [GeV]; entries/bin', 25, 0.0, 200.0)
        h_bestleta= r.TH1F('h_bestleta_'+sample, '; #eta, best soft lepton; entries/bin', 25, -3.0, +3.0)
        h_bestdphi= r.TH1F('h_bestdphi'+sample, '; |#Delta#phi(l_{hard},l_{soft})|, best pair [rad]; entries/bin',
                           25, 0.0, +2.*math.pi)
        h_bestdrj = r.TH1F('h_bestdrj'+sample, '; min#DeltaR(l_{soft}), best soft lep; entries/bin',
                           20, 0.0, 2.0)
        for event in tree :
            l0 = addTlv(event.l0)
            l1 = addTlv(event.l1)
            met = addTlv(event.met)
            jets = [addTlv(j) for j in event.jets]
            lepts = [addTlv(l) for l in event.lepts]
            pars = event.pars
            weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
            h_nlep.Fill(float(len(lepts)), weight)
            if len(lepts) < 1 : continue
            lepts = filter(lambda l : lepIsSeparatedFromOther(l.p4, [l0.p4, l1.p4]), lepts)
            h_nillep.Fill(float(len(lepts)), weight)
            validPairs = [(lh, ls) for lh in [l0, l1] for ls in lepts if lepPairIsZcand(lh, ls)]
            h_npairs.Fill(float(len(validPairs)), weight)
            pairsSortedByBestM = sorted(validPairs, key=deltaMZ0)
            for i,(lh,ls) in enumerate(pairsSortedByBestM) :
                m = (lh.p4+ls.p4).M()
                h_masses.Fill(m, weight)
                if i==0 :
                    h_bestm.Fill(m, weight)
                    h_bestlpt.Fill(ls.p4.Pt(), weight)
                    h_bestleta.Fill(ls.p4.Eta(), weight)
                    h_bestdphi.Fill(abs(phi_mpi_pi(ls.p4.DeltaPhi(lh.p4))), weight)
                    if len(jets) :
                        minDrJ = sorted([ls.p4.DeltaR(j.p4) for j in jets])[0]
                        h_bestdrj.Fill(minDrJ, weight)
            if iEvent<10 : print "jets [%d], lepts[%d]"%(len(jets), len(lepts))
            iEvent += 1
        hs_nlep    [sample] = h_nlep
        hs_nillep  [sample] = h_nillep
        hs_npairs  [sample] = h_npairs
        hs_masses  [sample] = h_masses
        hs_bestm   [sample] = h_bestm
        hs_bestlpt [sample] = h_bestlpt
        hs_bestleta[sample] = h_bestleta
        hs_bestdphi[sample] = h_bestdphi
        hs_bestdrj [sample] = h_bestdrj
    return {'hs_nlep'    : hs_nlep,
            'hs_nillep'  : hs_nillep,
            'hs_npairs'  : hs_npairs,
            'hs_masses'  : hs_masses,
            'hs_bestm'   : hs_bestm,
            'hs_bestlpt' : hs_bestlpt,
            'hs_bestleta': hs_bestleta,
            'hs_bestdphi': hs_bestdphi,
            'hs_bestdrj' : hs_bestdrj
            }

def plot(histos, var) :
    c = r.TCanvas('c_'+var,'')
    c.cd()
    pm = first(histos)
    pm.SetMaximum(1.1*max([h.GetMaximum() for h in histos.values()]))
    pm.Draw('axis')
    markers = dict(zip(histos.keys(),
                       [r.kPlus, r.kCircle, r.kMultiply, r.kOpenSquare, r.kOpenTriangleUp, r.kOpenTriangleDown]))
    for s,h in histos.iteritems() :
        h.SetLineColor(colors[s] if s in colors else r.kBlack)
        h.SetMarkerColor(h.GetLineColor())
        h.SetLineWidth(2*h.GetLineWidth())
        h.SetMarkerStyle(markers[s])
        h.Draw('same')
    leg = drawLegendWithDictKeys(c, histos)
    leg.Draw()
    c.Update()
    c.SaveAs(var+'.png')


if __name__=='__main__' :
    histos = fillHistos()
    for v, hs in histos.iteritems() :
        plot(hs, v)
