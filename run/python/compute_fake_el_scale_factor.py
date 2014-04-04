#!/bin/env python


import array
import glob
import os
from utils import first
import rootUtils
from rootUtils import importRoot, importRootCorePackages
r = rootUtils.importRoot()
r.gROOT.SetStyle('Plain')
rootUtils.importRootCorePackages()
from SampleUtils import colors


def main():
    tag = 'Apr_02'
    inputDir = 'out/fakerate'
    samples = ['Pythia8B_AU2_CTEQ6L1_bbTomu20',
               'Pythia8B_AU2_CTEQ6L1_ccTomu20'
               ]
    inputFileNames = [inputDir+"/%(sample)s_fake_tuple_%(tag)s.root"%{'sample':s, 'tag':tag} for s in samples]
    treeName = 'HeavyFlavorControlRegion'
    chain = r.TChain(treeName)
    [chain.Add(ifn) for ifn in inputFileNames]
    print treeName
    print '\n'.join(inputFileNames)
    print "about to loop on %d entries"%chain.GetEntries()
#     h_slpt = r.TH1F('h_slpt_'+sample, '; p_{T, soft-lep}[GeV]; entries/bin', 35, 0.0, 700.0)
    totNelec, totWeightLoose, totWeightTight = 0, 0.0, 0.0
    nElecTight = 0
    nOutRange, wOutRange = 0, 0.0
    for event in chain :
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        probe = event.l1
        isEl, isTight = probe.isEl, probe.isTight
        tlv = r.TLorentzVector()
        tlv.SetPxPyPzE(probe.px, probe.py, probe.pz, probe.E)
        pt = tlv.Pt()
        if not isEl : continue
        if pt<10.0 or pt>100.0 :
            nOutRange += 1
            wOutRange += weight
        totNelec += 1
        totWeightLoose += weight
        if isTight :
            totWeightTight += weight
            nElecTight += 1
            print 'tight electron ',runN,' ',evtN,' px ',probe.px
        
    print 'totNelec : ',totNelec
    print 'nElecTight : ',nElecTight
    print 'totWeightLoose: ',totWeightLoose
    print 'totWeightTight: ',totWeightTight
    print 'nOutRange : ',nOutRange
    print 'wOutRange : ',wOutRange
    histoNames = ['elec_fakeHF_high_all_l_pt_num',
                  'elec_fakeHF_high_all_l_pt_den',
                  ]
    for sample, tupleFname in zip(samples, inputFileNames) :
        histoFname = tupleFname.replace('_fake_tuple','')
        print histoFname
        histoFile = r.TFile.Open(histoFname)
        histos = [histoFile.Get(hn) for hn in histoNames]
        print sample+' : '+' '.join(["%s : %.2f (%d)"%(h.GetName(), h.Integral(0, -1), h.GetEntries()) for h in histos])

if __name__=='__main__':
    main()
