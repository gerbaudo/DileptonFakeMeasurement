#!/bin/env python

# Generic look at the truth jets for the higgsino signal samples
#
# davide.gerbaudo@gmail.com
# Feb 2013

import itertools, os, sys
import ROOT as r
import gen

def containerDirectory():
    "Full path of the directory where this file is located"
    path = os.path.realpath(__file__)
    return os.path.dirname(path)

r.gROOT.SetBatch(1)
r.gStyle.SetOptStat('nemri')
r.gSystem.Load(containerDirectory()+'/../../../SusyNtuple/StandAlone/libSusyNtuple.so')
r.gInterpreter.GenerateDictionary('vector< vector< int > >', 'vector')

treeName = 'susy'
inputFileName = '/tmp/gerbaudo/'\
                'mc12_8TeV.176575.Herwigpp_UEEE3_CTEQ6L1_simplifiedModel_wA_noslep_WH_2Lep_2.merge.NTUP_SUSY.e1702_s1581_s1586_r3658_r3549_p1328_tid01176849_00/'\
                'NTUP_SUSY.01176849._000001.root.1'

if len(sys.argv) > 1 : inputFileName = sys.argv[1]

inputFile = r.TFile.Open(inputFileName)
tree = inputFile.Get(treeName)
nEntries = tree.GetEntries()
#nEntries = 1000

GeV = 1.0
MeV2GeV = 0.001
tlv = r.TLorentzVector

minJetPt = 10.0*GeV
hTruthNjets = r.TH1F('hTruthNjets', 'truth jets (p_{T}>%dGeV); N_{jets}; jets/event'%minJetPt, 11, -0.5, 10.5)
hTruthJetEta = r.TH1F('hTruthJetEta', 'truth jets #eta; #eta; jets', 25, 0., 5.)
hTruthJetPt  = r.TH1F('hTruthJetPt',  'truth jets p_{T}; p_{T}; jets', 25, 0., 250.)
hTruthJet0Pt  = r.TH1F('hTruthJet0Pt',  'truth 1st jet p_{T}; p_{T}; jets', 25, 0., 250.)
hTruthJet1Pt  = r.TH1F('hTruthJet1Pt',  'truth 2nd jet p_{T}; p_{T}; jets', 25, 0., 250.)
hTruthJet0Eta  = r.TH1F('hTruthJet0Eta',  'truth 1st jet #eta; #eta; jets', 25, 0., 5.)
hTruthJet1Eta  = r.TH1F('hTruthJet1Eta',  'truth 2nd jet #eta; #eta; jets', 25, 0., 5.)
ptThresholds = [10., 15, 20., 25, 30.]
ptColors = [r.kRed, r.kBlue, r.kOrange, r.kMagenta, r.kAzure, ]
ptLineStyles = int(len(ptColors)/2+1)*[1, 2]
ptLineWidths = int(len(ptColors)/2+1)*[2, 3]
hsNjetAbovePt  = [r.TH1F('hNjetAbove%d'%p,  'N jets with p_{T}>%d; N_{jets}; jets/event'%p, 11, -0.5, 10.5) for p in ptThresholds]
hsMassJetPair  = [r.TH1F('hMassJetPair_%d_%d_%dj'%(i,j,nJet),  'dijet mass (%d, %d), %d jets'%(i,j, nJet), 25, 0.0, 250.)
                  for nJet in [2,3] # consider only the cases with 2 or 3 jets
                  for i,j in itertools.combinations(range(nJet), 2)]
hDecays = ['WW', 'ZZ', 'tautau', 'bbar', 'mumu', 'unknown']
hsNjetAbovePtPerHdecay = dict(zip(hDecays,
                                  [[r.TH1F('hNjetAbove%d_%s'%(p,d), 'N jets with p_{T}>%d, H #to %s; N_{jets}; jets/event'%(p,d), 11, -0.5, 10.5)
                                    for p in ptThresholds]
                                   for d in hDecays]))

print "looping over %d entries"%nEntries
findHiggs = gen.findInterestingHiggsWithChiAndPar
for iEntry in xrange(nEntries) :
    tree.GetEntry(iEntry)
    nTruthJets = tree.jet_AntiKt4TruthJets_n
    truthJets = [tlv(0., 0., 0., 0.) for i in xrange(nTruthJets)]
    pdg, parents, children = tree.mc_pdgId, tree.mc_parent_index, tree.mc_child_index
    interestingIhiggs, higgsChildren, higgsParents = findHiggs(pdg, parents, children)
    if len(interestingIhiggs)>1 :
        print "skip event with multiple (%d) higgs"%len(interestingIhiggs)
        continue
    higgsChildren, higgsParents = higgsChildren[0], higgsParents[0]
    hDecay = gen.guessHdecayLabel(higgsChildren)
    for j, pt, eta, phi, en in zip(truthJets,
                                   tree.jet_AntiKt4TruthJets_pt,
                                   tree.jet_AntiKt4TruthJets_eta,
                                   tree.jet_AntiKt4TruthJets_phi,
                                   tree.jet_AntiKt4TruthJets_E) :
        j.SetPtEtaPhiE(pt*MeV2GeV, eta, phi, en*MeV2GeV)
    truthJets = [j for j in truthJets if j.Pt()>minJetPt]
    nTruthJets = len(truthJets)
    sorted(truthJets, key=lambda j : j.Pt())
    hTruthNjets.Fill(len(truthJets))
    for i, j in enumerate(truthJets) :
        hTruthJetEta.Fill(j.Eta())
        hTruthJetPt.Fill(j.Pt())
        if   i==0 :
            hTruthJet0Pt.Fill(j.Pt())
            hTruthJet0Eta.Fill(j.Eta())
        elif i==1 :
            hTruthJet1Pt.Fill(j.Pt())
            hTruthJet1Eta.Fill(j.Eta())
    for i,t,h in zip(range(len(ptThresholds)), ptThresholds, hsNjetAbovePt) :
        nj = len([j for j in truthJets if j.Pt()>t])
        h.Fill(nj)
        hsNjetAbovePtPerHdecay[hDecay][i].Fill(nj)
    nJetsForPairs = 3 if nTruthJets>2 else (2 if nTruthJets==2 else 0)
    for i, (iJ,jJ) in enumerate(itertools.combinations(range(nJetsForPairs), 2)) :
        m_ij = (truthJets[iJ] + truthJets[jJ]).M()
        hsMassJetPair[(0 if nJetsForPairs==2 else 1) + i].Fill(m_ij)        
        
for h in [hTruthNjets, hTruthJetEta, hTruthJetPt, hTruthJet0Pt, hTruthJet1Pt, hTruthJet0Eta, hTruthJet1Eta] :
    print "%s(%d) , mean %.1f, RMS %.1f"%(h.GetName(), h.GetEntries(), h.GetMean(), h.GetRMS())
    #print [h.GetBinContent(i) for i in range(0,h.GetNbinsX()+2)]
    can = r.TCanvas('c_'+h.GetName(),'')
    can.cd()
    h.Draw()
    can.SaveAs(can.GetName()+'.png')

for h in hsMassJetPair :
    print "%s(%d) , mean %.1f, RMS %.1f"%(h.GetName(), h.GetEntries(), h.GetMean(), h.GetRMS())
    can = r.TCanvas('c_'+h.GetName(),'')
    can.cd()
    h.Draw()
    can.SaveAs(can.GetName()+'.png')

maxEtaCut = 2.5
for h in [hTruthJetEta, hTruthJet0Eta, hTruthJet1Eta] :
    totInt = h.Integral()
    partInt = h.Integral(h.FindBin(maxEtaCut), h.GetNbinsX())
    if totInt : print "%s : above %.1f -> %.2f"%(h.GetName(), maxEtaCut, 100.*partInt/totInt)



def drawPtThresHistos(histos, label,
                      ptThresholds=ptThresholds, lColors=ptColors, lStyles=ptLineStyles, lWidths=ptLineWidths) :
    can = r.TCanvas('c_hsNjetAbovePt'+label, '')
    can.cd()
    padMaster = histos[0].Clone('padMaster'+label)
    padMaster.Clear()
    padMaster.SetMaximum(1.1*max([h.GetMaximum() for h in histos]))
    padMaster.SetMinimum(1.0*min([h.GetMinimum() for h in histos]))
    padMaster.SetStats(0)
    padMaster.Draw()
    leg = r.TLegend(0.65, 0.65, 0.9, 0.9, label)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    for h, c, s, w in zip(histos, lColors, lStyles, lWidths) :
        h.SetLineColor(c)
        h.SetLineStyle(s)
        h.SetLineWidth(w)
    for h,t in zip(histos, ptThresholds) :
        h.Draw('same')
        leg.AddEntry(h,"p_{T}>%dGeV, <N>=%.2f"%(t, h.GetMean()), 'l')
    leg.Draw()
    can.SaveAs(can.GetName()+'.png')
for label, histos in [('any', hsNjetAbovePt)] + [(l,h) for l,h in hsNjetAbovePtPerHdecay.iteritems()]:
    drawPtThresHistos(histos, label)

# list of truth-jet branches, just for ref
#
# jet_AntiKt4TruthJets_n : Int_t
# jet_AntiKt4TruthJets_E : vector<float>
# jet_AntiKt4TruthJets_pt : vector<float>                            *
# jet_AntiKt4TruthJets_m : vector<float>                             *
# jet_AntiKt4TruthJets_eta : vector<float>                           *
# jet_AntiKt4TruthJets_phi : vector<float>                           *
# jet_AntiKt4TruthJets_EtaOrigin : vector<float>                     *
# jet_AntiKt4TruthJets_PhiOrigin : vector<float>                     *
# jet_AntiKt4TruthJets_MOrigin : vector<float>                       *
# jet_AntiKt4TruthJets_WIDTH : vector<float>                         *
# jet_AntiKt4TruthJets_n90 : vector<float>                           *
# jet_AntiKt4TruthJets_Timing : vector<float>                        *
# jet_AntiKt4TruthJets_LArQuality : vector<float>                    *
# jet_AntiKt4TruthJets_nTrk : vector<float>                          *
# jet_AntiKt4TruthJets_sumPtTrk : vector<float>                      *
# jet_AntiKt4TruthJets_OriginIndex : vector<float>                   *
# jet_AntiKt4TruthJets_HECQuality : vector<float>                    *
# jet_AntiKt4TruthJets_NegativeE : vector<float>                     *
# jet_AntiKt4TruthJets_AverageLArQF : vector<float>                  *
# jet_AntiKt4TruthJets_BCH_CORR_CELL : vector<float>                 *
# jet_AntiKt4TruthJets_BCH_CORR_DOTX : vector<float>                 *
# jet_AntiKt4TruthJets_BCH_CORR_JET : vector<float>                  *
# jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL : vector<float>          *
# jet_AntiKt4TruthJets_ENG_BAD_CELLS : vector<float>                 *
# jet_AntiKt4TruthJets_N_BAD_CELLS : vector<float>                   *
# jet_AntiKt4TruthJets_N_BAD_CELLS_CORR : vector<float>              *
# jet_AntiKt4TruthJets_BAD_CELLS_CORR_E : vector<float>              *
# jet_AntiKt4TruthJets_NumTowers : vector<float>                     *
# jet_AntiKt4TruthJets_ootFracCells5 : vector<float>                 *
# jet_AntiKt4TruthJets_ootFracCells10 : vector<float>                *
# jet_AntiKt4TruthJets_ootFracClusters5 : vector<float>              *
# jet_AntiKt4TruthJets_ootFracClusters10 : vector<float>             *
# jet_AntiKt4TruthJets_SamplingMax : vector<int>                     *
# jet_AntiKt4TruthJets_fracSamplingMax : vector<float>               *
# jet_AntiKt4TruthJets_hecf : vector<float>                          *
# jet_AntiKt4TruthJets_tgap3f : vector<float>                        *
# jet_AntiKt4TruthJets_isUgly : vector<int>                          *
# jet_AntiKt4TruthJets_isBadLooseMinus : vector<int>                 *
# jet_AntiKt4TruthJets_isBadLoose : vector<int>                      *
# jet_AntiKt4TruthJets_isBadMedium : vector<int>                     *
# jet_AntiKt4TruthJets_isBadTight : vector<int>                      *
# jet_AntiKt4TruthJets_emfrac : vector<float>                        *
# jet_AntiKt4TruthJets_Offset : vector<float>                        *
# jet_AntiKt4TruthJets_EMJES : vector<float>                         *
# jet_AntiKt4TruthJets_EMJES_EtaCorr : vector<float>                 *
# jet_AntiKt4TruthJets_EMJESnooffset : vector<float>                 *
# jet_AntiKt4TruthJets_LCJES : vector<float>                         *
# jet_AntiKt4TruthJets_LCJES_EtaCorr : vector<float>                 *
# jet_AntiKt4TruthJets_emscale_E : vector<float>                     *
# jet_AntiKt4TruthJets_emscale_pt : vector<float>                    *
# jet_AntiKt4TruthJets_emscale_m : vector<float>                     *
# jet_AntiKt4TruthJets_emscale_eta : vector<float>                   *
# jet_AntiKt4TruthJets_ActiveArea : vector<float>                    *
# jet_AntiKt4TruthJets_ActiveAreaPx : vector<float>                  *
# jet_AntiKt4TruthJets_ActiveAreaPy : vector<float>                  *
# jet_AntiKt4TruthJets_ActiveAreaPz : vector<float>                  *
# jet_AntiKt4TruthJets_ActiveAreaE : vector<float>                   *
# jet_AntiKt4TruthJets_flavor_truth_label : vector<int>              *
# jet_AntiKt4TruthJets_flavor_truth_dRminToB : vector<float>         *
# jet_AntiKt4TruthJets_flavor_truth_dRminToC : vector<float>         *
# jet_AntiKt4TruthJets_flavor_truth_dRminToT : vector<float>         *
# jet_AntiKt4TruthJets_flavor_truth_BHadronpdg : vector<int>         *
# jet_AntiKt4TruthJets_flavor_truth_vx_x : vector<float>             *
# jet_AntiKt4TruthJets_flavor_truth_vx_y : vector<float>             *
# jet_AntiKt4TruthJets_flavor_truth_vx_z : vector<float>             *
# jet_AntiKt4TruthJets_el_dr : vector<float>                         *
# jet_AntiKt4TruthJets_el_matched : vector<int>                      *
# jet_AntiKt4TruthJets_mu_dr : vector<float>                         *
# jet_AntiKt4TruthJets_mu_matched : vector<int>                      *
