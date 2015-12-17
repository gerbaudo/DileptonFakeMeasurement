#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <sstream>      // std::ostringstream
#include <stdexcept>
#include "TCanvas.h"
#include "DileptonFakeMeasurement/SusySelection.h"

#include "DileptonFakeMeasurement/criteria.h"
#include "DileptonFakeMeasurement/kinematic.h"
#include "DileptonFakeMeasurement/utils.h"


using namespace std;
using namespace Susy;
namespace swk = susy::wh::kin;

std::string SusySelection::WeightComponents::str() const
{
  std::ostringstream oss;
  oss<<" susynt: "<<susynt
     <<" lepSf: "<<lepSf
     <<" btag: "<<btag
     <<" trigger: "<<trigger
     <<" qflip: "<<qflip
     <<" fake: "<<fake;
  return oss.str();
}
//-----------------------------------------
SusySelection::SusySelection() :
  m_tupleMaker("",""),
  m_writeTuple(false),
  m_debugThisEvent(false),
  m_outTupleFile(""),
  // m_trigObj(NULL),
  m_useMCTrig(false),
  m_w(1.0),
  m_qflipProb(0.0)
{
  resetAllCounters();
//  setAnaType(Ana_2LepWH); // tmp test for HLFV
  setAnaType(Ana_2Lep);
  setSelectTaus(true);
  initChargeFlipTool();
}
void SusySelection::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  if(m_dbg) cout << "SusySelection::Begin" << endl;
  string period = "Moriond";
  bool useReweightUtils = false;
  // m_trigObj = new DilTrigLogic(period, useReweightUtils);
  // if(m_useMCTrig) m_trigObj->useMCTrigger();
  if(m_writeTuple) {
      if(endswith(m_outTupleFile, ".root") && m_tupleMaker.init(m_outTupleFile, "SusySel"))
          cout<<"initialized ntuple file "<<m_outTupleFile<<endl;
      else {
          cout<<"cannot initialize ntuple file '"<<m_outTupleFile<<"'"<<endl;
          m_writeTuple = false;
      }
  }
}
//-----------------------------------------
JetVector SusySelection::filterClJets(const JetVector &jets)
{
    JetVector oj;
    for(size_t i=0; i<jets.size(); ++i)
        if(SusyNtTools::isCentralLightJet(jets[i])) oj.push_back(jets[i]);
    return oj;
}
//-----------------------------------------
Bool_t SusySelection::Process(Long64_t entry)
{
  m_printer.countAndPrint(cout);
  GetEntry(entry);
  clearObjects();
  cacheStaticWeightComponents();
  increment(n_readin, m_weightComponents);
  bool removeLepsFromIso(false), allowQflip(true);
  selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium);
  if(!selectEvent()) return kTRUE;
  const JetVector&   bj = m_baseJets;
  const LeptonVector& l = m_signalLeptons;
  if(l.size()>1) computeNonStaticWeightComponents(l, bj); else return false;
  if(passSrSs(WH_SRSS1,
              m_signalLeptons, m_signalTaus, m_signalJets2Lep, m_met, allowQflip)
     .lepPt) {
      if(m_writeTuple) {
          double weight(m_weightComponents.product());
          unsigned int run(nt.evt()->run), event(nt.evt()->event);
          LeptonVector anyLep(getAnyElOrMu(nt));
          LeptonVector lowPtLep(subtract_vector(anyLep, m_baseLeptons));
          const Lepton *l0 = m_signalLeptons[0];
          const Lepton *l1 = m_signalLeptons[1];
          const JetVector clJets(filterClJets(m_signalJets2Lep));
          m_tupleMaker.fill(weight, run, event, *l0, *l1, *m_met, lowPtLep, clJets);
      }
  }
  checkAndIncrementEwk(m_signalLeptons, m_signalJets2Lep, m_met);
  return kTRUE;
}
//-----------------------------------------
void SusySelection::Terminate()
{
    if(m_writeTuple) m_tupleMaker.close();
  SusyNtAna::Terminate();
  if(m_dbg) cout << "SusySelection::Terminate" << endl;
  dumpEventCounters();
  if(m_chargeFlip) delete m_chargeFlip;
}
//-----------------------------------------
void SusySelection::increment(float counters[], const WeightComponents &wc)
{
  counters[kRaw ] += 1.0;
  counters[kEvt ] += wc.gen;
  counters[kPU  ] += wc.gen * wc.pileup;
  counters[kLSF ] += wc.gen * wc.lepSf;
  counters[kBtag] += wc.gen * wc.btag;
  counters[kTrig] += wc.gen * wc.trigger;
  counters[kAll ] += wc.product();
}
//-----------------------------------------
bool SusySelection::selectEvent()
{
  if(m_dbg) cout << "SusySelection::selectEvent" << endl;
  // Basic event cuts
  int flag = nt.evt()->cutFlags[NtSys_NOM];
  //int hdec = nt.evt()->hDecay;
  const LeptonVector& bleps = m_baseLeptons;
  const JetVector &jets = m_baseJets;
  const JetVector &pjets = m_preJets;
  const Susy::Met *met = m_met;
  uint run = nt.evt()->run;
  bool mc = nt.evt()->isMC;
  m_debugThisEvent = susy::isEventInList(nt.evt()->event);
  float mllMin(20);
  WeightComponents &wc = m_weightComponents;
  if(passGRL        (flag           ))  { increment(n_pass_Grl     , wc);} else { return false; }
  if(passLarErr     (flag           ))  { increment(n_pass_LarErr  , wc);} else { return false; }
  if(passTileErr    (flag           ))  { increment(n_pass_TileErr , wc);} else { return false; }
  if(passTTCVeto    (flag           ))  { increment(n_pass_TTCVeto , wc);} else { return false; }
  if(passGoodVtx    (flag           ))  { increment(n_pass_GoodVtx , wc);} else { return false; }
  if(passTileTripCut(flag           ))  { increment(n_pass_TileTrip, wc);} else { return false; }
  if(passLAr        (flag           ))  { increment(n_pass_LAr     , wc);} else { return false; }
  if(!hasBadJet     (jets           ))  { increment(n_pass_BadJet  , wc);} else { return false; }
  if(passDeadRegions(pjets,met,run,mc)) { increment(n_pass_FEBCut  , wc);} else { return false; }
  if(!hasBadMuon    (m_preMuons     ))  { increment(n_pass_BadMuon , wc);} else { return false; }
  if(!hasCosmicMuon (m_baseMuons    ))  { increment(n_pass_Cosmic  , wc);} else { return false; }
  if(passHfor       (               ))  { increment(n_pass_hfor    , wc);} else { return false; }
  //if(passHtautauVeto(hdec)) { increment(n_pass_HttVeto ); } else { return false; }
  if(bleps.size() >= 2               )  { increment(n_pass_ge2l    , wc);} else { return false; }
  if(bleps.size() == 2               )  { increment(n_pass_eq2l    , wc);} else { return false; }
  if(susy::passMllMin(bleps, mllMin ))  { increment(n_pass_mll     , wc);} else { return false; }
  return true;
}
//-----------------------------------------
bool SusySelection::passSR6base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return susy::oppositeSign(leptons) && susy::sameFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR7base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return susy::oppositeSign(leptons) && susy::oppositeFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR8base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return susy::sameSign(leptons) && susy::sameFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR9base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return susy::sameSign(leptons) && susy::oppositeFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR6(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return false; // now obsolete
}
//-----------------------------------------
bool SusySelection::passSR7(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return false; // now obsolete
}
//-----------------------------------------
bool SusySelection::passSR8(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return false; // now obsolete
}
//-----------------------------------------
bool SusySelection::passSR9(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return false; // now obsolete
}
//-----------------------------------------
bool SusySelection::passSrSsBase()
{
  return false;
}
//-----------------------------------------
SsPassFlags SusySelection::passSrSs(const WH_SR signalRegion,
                                    LeptonVector& leptons,
                                    const TauVector& taus,
                                    const JetVector& jets,
                                    const Met *met,
                                    bool allowQflip)
{
  SsPassFlags f;
  if(leptons.size()<2) return f;
  DiLepEvtType ll(getDiLepEvtType(leptons));
  if(ll==ET_me) ll = ET_em;
  bool update4mom(true); // charge flip
  bool mc(nt.evt()->isMC), data(!mc);
  const LeptonVector &ls = leptons;
  LeptonVector &ncls = leptons; // non-const leptons: can be modified by qflip
  Met ncmet(*m_met); // non-const met \todo: should modify a non-const input
  const JetVector    &js = jets;
  WeightComponents &wc = m_weightComponents;

  increment(n_pass_category[ll], wc);
  if(susy::passNlepMin(ls, 2))                  { increment(n_pass_nSigLep  [ll], wc); f.eq2l       =true;} else return f;
  if(m_signalTaus.size()==0)                    { increment(n_pass_tauVeto  [ll], wc); f.tauVeto    =true;} else return f;
  if(passTrig2L     (ls))                       { increment(n_pass_tr2L     [ll], wc); f.trig2l     =true;} else return f;
  if(passTrig2LMatch(ls))                       { increment(n_pass_tr2LMatch[ll], wc); f.trig2lmatch=true;} else return f;
  if(data || susy::isTrueDilepton(ls))          { increment(n_pass_mcTrue2l [ll], wc); f.true2l     =true;} else return f;
  bool sameSign = allowQflip ? sameSignOrQflip(ncls, ncmet, ll, update4mom, mc) : susy::sameSign(ncls);
  if(sameSign)                                  { increment(n_pass_ss       [ll], wc); f.sameSign   =true;} else return f;
  met = &ncmet; // after qflip, use potentially smeared lep and met
  increment(n_pass_muIso    [ll], wc);
  increment(n_pass_elD0Sig  [ll], wc);
  f = assignNjetFlags(js, f);
  LeptonVector anyLeptons(getAnyElOrMu(nt));
  LeptonVector lowPtLep(subtract_vector(anyLeptons, m_baseLeptons));
  /*const*/ swk::DilepVars v(swk::compute2lVars(leptons, met, jets));
  v.l3veto = f.veto3rdL = passThirdLeptonVeto(ncls[0], ncls[1], lowPtLep, m_debugThisEvent);
  if(f.veto3rdL) increment(n_pass_3rdLep [ll], wc); else return f;
  if(f.fjveto  ) increment(n_pass_fjVeto [ll], wc); else return f;
  if(f.bjveto  ) increment(n_pass_bjVeto [ll], wc); else return f;
  if(f.ge1j    ) increment(n_pass_ge1j   [ll], wc); else return f;
  if(f.eq1j    ) increment(n_pass_eq1j   [ll], wc);
  else           increment(n_pass_ge2j   [ll], wc);
  assert(f.eq1j != f.ge2j); // from here on count separately eq1j and ge2j
  if(f.eq1j) SusySelection::passSrWh1j(v, f);
  else       SusySelection::passSrWh2j(v, f);
  float *cnt_lepPt    = f.eq1j ? n_pass_eq1jlepPt    [ll] : n_pass_ge2jlepPt    [ll];
  float *cnt_mllZveto = f.eq1j ? n_pass_eq1jmllZveto [ll] : n_pass_ge2jmllZveto [ll];
  float *cnt_detall   = f.eq1j ? n_pass_eq1jDetall   [ll] : n_pass_ge2jDetall   [ll];
  float *cnt_mtmax    = f.eq1j ? n_pass_eq1jMtMax    [ll] : n_pass_ge2jMtMax    [ll];
  float *cnt_mljj     = f.eq1j ? n_pass_eq1jMlj      [ll] : n_pass_ge2jMljj     [ll];
  float *cnt_ht       = f.eq1j ? n_pass_eq1jht       [ll] : n_pass_ge2jht       [ll];
  float *cnt_metRel   = f.eq1j ? n_pass_eq1jmetRel   [ll] : n_pass_ge2jmetRel   [ll];
  float *cnt_mWwt     = f.eq1j ? n_pass_eq1jmWwt     [ll] : n_pass_ge2jmWwt     [ll];
  if(f.lepPt   ) increment(cnt_lepPt   , wc); else return f;
  if(f.zllVeto ) increment(cnt_mllZveto, wc); else return f;
  if(f.dEtall  ) increment(cnt_detall  , wc); else return f;
  if(f.maxMt   ) increment(cnt_mtmax   , wc); else return f;
  if(f.mljj    ) increment(cnt_mljj    , wc); else return f;
  if(f.ht      ) increment(cnt_ht      , wc); else return f;
  if(f.metrel  ) increment(cnt_metRel  , wc); else return f;
  if(f.mtllmet ) increment(cnt_mWwt    , wc); else return f;
  return f;
}
//-----------------------------------------
bool SusySelection::passHfor(Susy::SusyNtObject &nto)
{
    // DG : inheriting hardcoded magic values from HforToolD3PD.cxx, dah.
    const int kill(4);
    return nto.evt()->hfor != kill;
}
//-----------------------------------------
bool SusySelection::passTrig2L(const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt)
{
  if(leptons.size() != 2 || !dtl) return false;
  return dtl->passDilEvtTrig(leptons, met, evt);
}
//-----------------------------------------
bool SusySelection::passTrig2LMatch(const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt)
{
  if(leptons.size() != 2 || !dtl) return false;
  return dtl->passDilTrigMatch(leptons, met, evt);
}
//-----------------------------------------
bool SusySelection::passTrig2LwithMatch(const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt)
{
  return (passTrig2L(leptons, dtl, met, evt) && passTrig2LMatch(leptons, dtl, met, evt));
}
//-----------------------------------------
bool SusySelection::sameSignOrQflip(LeptonVector& leptons, Met &met,
                                    const DiLepEvtType eventType,
                                    bool update4mom, bool isMc)
{
  bool isSS(susy::sameSign(leptons));
  if(isSS) return true;
  if(!isMc) return isSS;
  bool isOS(!isSS);
  bool canBeQflip(isOS && (eventType==ET_ee || eventType==ET_em || eventType==ET_me));
  if (!canBeQflip){ return false; }
  if(canBeQflip){
    uint systematic=NtSys_NOM; // DG sys todo
    m_qflipProb = computeChargeFlipProb(leptons, met, systematic, update4mom);
    m_weightComponents.qflip = m_qflipProb;
  }
  return true;
}
//-----------------------------------------
bool SusySelection::passJetVeto(const JetVector& jets)
{
  // Require no light, b, or forward jets
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 + N_B20 + N_F30 == 0);
}
//-----------------------------------------
bool SusySelection::passbJetVeto(const JetVector& jets)
{
  // Reject if there is a b jet using 2L definition
  int N_B20 = numberOfCBJets(jets);
  return (N_B20 == 0);
}
//-----------------------------------------
bool SusySelection::passfJetVeto(const JetVector& jets)
{
  return (0 == numberOfFJets(jets));
}
//-----------------------------------------
bool SusySelection::passge1Jet(const JetVector& jets)
{
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 >=1 && N_B20 + N_F30 == 0);
}
bool SusySelection::passge2Jet(const JetVector& jets)
{
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 >=2 && N_B20 + N_F30 == 0);
}
bool SusySelection::passge3Jet(const JetVector& jets)
{
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 >=3 && N_B20 + N_F30 == 0);
}
//-----------------------------------------
bool SusySelection::passeq2Jet(const JetVector& jets)
{
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 == 2 && N_B20 + N_F30 == 0);
}
//-----------------------------------------
bool SusySelection::passge2JetWoutFwVeto(const JetVector& jets)
{
  return (numberOfCLJets(jets) >= 2 && numberOfCBJets(jets) < 1);
}
//-----------------------------------------
bool SusySelection::passeq2JetWoutFwVeto(const JetVector& jets)
{
  return (numberOfCLJets(jets) == 2 && numberOfCBJets(jets) < 1);
}
//-----------------------------------------
bool SusySelection::passMetRelMin(const Met *met, const LeptonVector& leptons,
                                  const JetVector& jets, float minVal){
  float metrel = getMetRel(met,leptons,jets);
  return (minVal < metrel);
}
//----------------------------------------------------------
bool SusySelection::passNj(const JetVector& jets, int minNj, int maxNj)
{
  int nj(numberOfCLJets(jets));
  return (minNj < nj && nj <= maxNj
	  && numberOfCBJets(jets) < 1);
}
//-----------------------------------------
bool SusySelection::passMuonRelIso(const LeptonVector &leptons, float maxVal)
{
  for(size_t i=0; i<leptons.size(); ++i){
    const Susy::Lepton* l = leptons[i];
    if(l->isMu()){
      const Muon* mu = static_cast<const Muon*>(l);
      if(!mu) continue;
      float etcone30 = muEtConeCorr(mu, m_baseElectrons, m_baseMuons,
                                    nt.evt()->nVtx, nt.evt()->isMC);
      if(mu->Pt() && (etcone30/mu->Pt() > maxVal)) return false;
    } // end if(isMu)
  } // end for(i)
  return true;
}
//-----------------------------------------
bool SusySelection::passEwkSs(const LeptonVector& leptons, const JetVector& jets, const Met* met)
{
    if(leptons.size()<2) return false;
    bool noBjets(numberOfCBJets(jets)==0), noFwJets(numberOfFJets(jets)==0);
    bool someCentralJets(numberOfCLJets(jets)>0);
    const Lepton &l0 = *leptons[0], &l1 = *leptons[1];
    TLorentzVector ll(l0+l1);
    return (noBjets && noFwJets && someCentralJets
            && (getMetRel(met, leptons, jets)>50.0)
            && susy::sameSign(leptons)
            && (ll.M()<60.0) && (ll.Pt()<20.) && (fabs(l0.DeltaPhi(l1)) >= 1.3));
}
//-----------------------------------------
bool SusySelection::passEwkSsLoose(const LeptonVector& leptons, const JetVector& jets, const Met* met)
{
    if(leptons.size()<2) return false;
    bool noBjets(numberOfCBJets(jets)==0), noFwJets(numberOfFJets(jets)==0);
    bool someCentralJets(numberOfCLJets(jets)>0);
    bool isEe(leptons[0]->isEle() && leptons[1]->isEle());
    bool passZeeVeto(isEe ? susy::passZllVeto(leptons, 91.2-10.0, 91.2+10.0) : true);
    return (noBjets && noFwJets && someCentralJets && passZeeVeto
            && susy::sameSign(leptons)
            && (getMetRel(met, leptons, jets)>40.0));
}
//-----------------------------------------
bool SusySelection::passEwkSsLea(const LeptonVector& leptons, const JetVector& jets, const Met* met)
{
    bool pass=false;
    if(leptons.size()>=2) {
        bool noBjets(numberOfCBJets(jets)==0), noFwJets(numberOfFJets(jets)==0);
        bool someCentralJets(numberOfCLJets(jets)>0);
        const Susy::Lepton &l0 = *leptons[0], &l1 = *leptons[1];
        float mll = (l0+l1).M();
        bool isEe(l0.isEle() && l1.isEle());
        bool passZeeVeto(isEe ? fabs(mll- 91.2)>10.0 : true);
        pass= (noBjets && noFwJets && someCentralJets && passZeeVeto
               && susy::sameSign(leptons)
               && (met->Et<40.0));
    }
    return pass;
}
//-----------------------------------------
void SusySelection::checkAndIncrementEwk(const LeptonVector& leptons, const JetVector& jets, const Met* met)
{
  DiLepEvtType ll(getDiLepEvtType(leptons));
  if(ll==ET_me) ll = ET_em;
  if(SusySelection::passEwkSs     (leptons, jets, met)) increment(n_pass_ewkSs      [ll], m_weightComponents);
  if(SusySelection::passEwkSsLoose(leptons, jets, met)) increment(n_pass_ewkSsLoose [ll], m_weightComponents);
  if(SusySelection::passEwkSsLea  (leptons, jets, met)) increment(n_pass_ewkSsLea   [ll], m_weightComponents);
}
//-----------------------------------------
void SusySelection::cacheStaticWeightComponents()
{
  m_weightComponents.reset();
  if(!nt.evt()->isMC) {m_weightComponents.reset(); return;}
  m_weightComponents.gen = nt.evt()->w;
  m_weightComponents.pileup = nt.evt()->wPileup;
  bool useSumwMap(true);
  m_weightComponents.susynt = (m_useXsReader ?
                               computeEventWeightXsFromReader(LUMI_A_L) :
                               SusyNtAna::getEventWeight(LUMI_A_L, useSumwMap));
  float genpu(m_weightComponents.gen*m_weightComponents.pileup);
  m_weightComponents.norm = (genpu != 0.0 ? m_weightComponents.susynt/genpu : 1.0);
}
//-----------------------------------------
void SusySelection::computeNonStaticWeightComponents(cvl_t& leptons, cvj_t& jets)
{
  if(!nt.evt()->isMC) {m_weightComponents.reset(); return;}
  m_weightComponents.lepSf = susy::getLeptonEff2Lep(leptons);
  m_weightComponents.trigger = getTriggerWeight2Lep(leptons);
  m_weightComponents.btag = getBTagWeight(jets, nt.evt());
}
//-----------------------------------------
float SusySelection::getBTagWeight(cvj_t& jets, const Event* evt)
{
  JetVector tempJets;
  for(uint ij=0; ij<jets.size(); ++ij){
    Jet* jet = jets.at(ij);
    if( !(jet->Pt() > 20 && fabs(jet->Eta()) < JET_ETA_CUT_2L) ) continue;
    tempJets.push_back(jet);
  }
  return bTagSF(evt, tempJets, evt->mcChannel, BTag_NOM);
}
//-----------------------------------------
float SusySelection::getTriggerWeight2Lep(const LeptonVector &leptons)
{
  float trigW = 1.0;
  return trigW;
  // // if m_useMCTrig, then we are dropping evts with DilTrigLogic::passDil*, not weighting them
  // // DG Jun2013 : verify this with Matt & Josephine
  // if(!m_useMCTrig){
  //   if(leptons.size()==2) trigW = m_trigObj->getTriggerWeight(leptons,
  //                                                             nt.evt()->isMC,
  //                                                             m_met->Et,
  //                                                             m_signalJets2Lep.size(),
  //                                                             nt.evt()->nVtx,
  //                                                             NtSys_NOM);
  //   bool twIsInvalid(isnan(trigW) || trigW<0.0);
  //   assert(!twIsInvalid);
  //   if(twIsInvalid){
  //     if(m_dbg)
  //       cout<<"SusySelection::getTriggerWeight: invalid weigth "<<trigW<<", using 0.0"<<endl;
  //     trigW = (twIsInvalid ? 0.0 : trigW);
  //   }
  // }
  // return trigW;
}
//-----------------------------------------
// helper function: write header with event types
std::string lineLabelsPerEventType(const string *labels, int colWidth){
  std::ostringstream oss;
  for(int i=0; i<ET_N-1; ++i)
    oss<<std::setw(colWidth)<<labels[i];
  oss<<std::setw(colWidth)<<"em+me";
  return oss.str();
}
// helper function: for a given weight type, write line with counts for each event type
std::string lineCountersPerEventType(const float cnt[ET_N][kWeightTypesN],
                                     int weightType, int colWidth){
  std::ostringstream oss;
  // bool raw(weightType==WT_Raw);
  // int precision(raw ? 0 : 2); // DG Aug2013 not working properly tobefixed
  for(int i=0; i<ET_N-1; ++i)
    oss<<std::setw(colWidth)
      //<<(raw ? std::fixed : "")
      //<<std::setprecision(precision)
       <<cnt[i][weightType];
  oss<<std::setw(colWidth)<<(cnt[ET_em][weightType] + cnt[ET_me][weightType]);
  return oss.str();
}
void SusySelection::dumpEventCounters()
{
  string v_ET[] = {"ee","mm","em","me"};
  string v_WT[] = {"Raw","Event","Pileup","LeptonSF","btagSF","TrigSF","All"};
  int colWidth(10);
  int &cw = colWidth;
  using std::setw;
  int nCols(ET_N-1);
  string topRule(nCols*colWidth, '*');
  string midRule(nCols*colWidth, '-');
  // define a function reference to shorten lines
  string (&lcpet)(const float cnt[ET_N][kWeightTypesN], int weightType, int colWidth) = lineCountersPerEventType;
  for(int w=0; w<kWeightTypesN; ++w){
    cout<<topRule                                                    <<endl;
    cout<<"Event counts for weight: "<< v_WT             [w]         <<endl;
    cout<<midRule                                                    <<endl;
    cout<<"input:           : "<<setw(cw)<<n_readin           [w]    <<endl;
    cout<<"GRL              : "<<setw(cw)<<n_pass_Grl         [w]    <<endl;
    cout<<"LarErr           : "<<setw(cw)<<n_pass_LarErr      [w]    <<endl;
    cout<<"TileErr          : "<<setw(cw)<<n_pass_TileErr     [w]    <<endl;
    cout<<"TTCVeto          : "<<setw(cw)<<n_pass_TTCVeto     [w]    <<endl;
    cout<<"GoodVtx          : "<<setw(cw)<<n_pass_GoodVtx     [w]    <<endl;
    cout<<"TileTripCut      : "<<setw(cw)<<n_pass_TileTrip    [w]    <<endl;
    cout<<"LAr:             : "<<setw(cw)<<n_pass_LAr         [w]    <<endl;
    cout<<"BadJet:          : "<<setw(cw)<<n_pass_BadJet      [w]    <<endl;
    cout<<"FEB:             : "<<setw(cw)<<n_pass_FEBCut      [w]    <<endl;
    cout<<"BadMu:           : "<<setw(cw)<<n_pass_BadMuon     [w]    <<endl;
    cout<<"Cosmic:          : "<<setw(cw)<<n_pass_Cosmic      [w]    <<endl;
    cout<<"hfor:            : "<<setw(cw)<<n_pass_hfor        [w]    <<endl;
    cout<<"Htautau veto     : "<<setw(cw)<<n_pass_HttVeto     [w]    <<endl;
    cout<<"atleast 2        : "<<setw(cw)<<n_pass_ge2l        [w]    <<endl;
    cout<<"exactly 2        : "<<setw(cw)<<n_pass_eq2l        [w]    <<endl;
    cout<<"mll              : "<<setw(cw)<<n_pass_mll         [w]    <<endl;
    cout<<"nSigLep          : "<<setw(cw)<<n_pass_signalLep   [w]    <<endl;
    cout<<"   ------  Start Comparison Here ------ "                 <<endl;
    cout<<"Dilepton type    : "<<lineLabelsPerEventType(v_ET, cw)    <<endl;
    cout<<"category         : "<<lcpet(n_pass_category       , w, cw)<<endl;
    cout<<"nSigLep          : "<<lcpet(n_pass_nSigLep        , w, cw)<<endl;
    cout<<"tauVeto          : "<<lcpet(n_pass_tauVeto        , w, cw)<<endl;
    cout<<"trig:            : "<<lcpet(n_pass_tr2L           , w, cw)<<endl;
    cout<<"trig match:      : "<<lcpet(n_pass_tr2LMatch      , w, cw)<<endl;
    cout<<"mc prompt2l      : "<<lcpet(n_pass_mcTrue2l       , w, cw)<<endl;
    cout<<"SS:              : "<<lcpet(n_pass_ss             , w, cw)<<endl;
    cout<<"muIso            : "<<lcpet(n_pass_muIso          , w, cw)<<endl;
    cout<<"elD0Sig          : "<<lcpet(n_pass_elD0Sig        , w, cw)<<endl;
    cout<<"3rdLepVeto       : "<<lcpet(n_pass_3rdLep         , w, cw)<<endl;
    cout<<"fjVeto           : "<<lcpet(n_pass_fjVeto         , w, cw)<<endl;
    cout<<"bjVeto           : "<<lcpet(n_pass_bjVeto         , w, cw)<<endl;
    cout<<"ge1j             : "<<lcpet(n_pass_ge1j           , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
    cout<<"eq1j             : "<<lcpet(n_pass_eq1j           , w, cw)<<endl;
    cout<<"lepPt            : "<<lcpet(n_pass_eq1jlepPt      , w, cw)<<endl;
    cout<<"mllZveto         : "<<lcpet(n_pass_eq1jmllZveto   , w, cw)<<endl;
    cout<<"Detall           : "<<lcpet(n_pass_eq1jDetall     , w, cw)<<endl;
    cout<<"MtMax            : "<<lcpet(n_pass_eq1jMtMax      , w, cw)<<endl;
    cout<<"Mlj              : "<<lcpet(n_pass_eq1jMlj        , w, cw)<<endl;
    cout<<"ht               : "<<lcpet(n_pass_eq1jht         , w, cw)<<endl;
    cout<<"metRel           : "<<lcpet(n_pass_eq1jmetRel     , w, cw)<<endl;
    cout<<"mWwt             : "<<lcpet(n_pass_eq1jmWwt       , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
    cout<<"ge2j             : "<<lcpet(n_pass_ge2j           , w, cw)<<endl;
    cout<<"lepPt            : "<<lcpet(n_pass_ge2jlepPt      , w, cw)<<endl;
    cout<<"mllZveto         : "<<lcpet(n_pass_ge2jmllZveto   , w, cw)<<endl;
    cout<<"Detall           : "<<lcpet(n_pass_ge2jDetall     , w, cw)<<endl;
    cout<<"MtMax            : "<<lcpet(n_pass_ge2jMtMax      , w, cw)<<endl;
    cout<<"Mljj             : "<<lcpet(n_pass_ge2jMljj       , w, cw)<<endl;
    cout<<"ht               : "<<lcpet(n_pass_ge2jht         , w, cw)<<endl;
    cout<<"metRel           : "<<lcpet(n_pass_ge2jmetRel     , w, cw)<<endl;
    cout<<"mWwt             : "<<lcpet(n_pass_ge2jmWwt       , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
    cout<<"   ------  Control regions ------ "                       <<endl;
    cout<<"cr1jzv           : "<<lcpet(n_pass_cr1jzv         , w, cw)<<endl;
    cout<<"cr1jfake         : "<<lcpet(n_pass_cr1jfake       , w, cw)<<endl;
    cout<<"cr1jzvfake       : "<<lcpet(n_pass_cr1jzvfake     , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
    cout<<"cr2jzv           : "<<lcpet(n_pass_cr2jzv         , w, cw)<<endl;
    cout<<"cr2jfake         : "<<lcpet(n_pass_cr2jfake       , w, cw)<<endl;
    cout<<"cr2jzvfake       : "<<lcpet(n_pass_cr2jzvfake     , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
    cout<<"ewkSs            : "<<lcpet(n_pass_ewkSs          , w, cw)<<endl;
    cout<<"ewkSsLoose       : "<<lcpet(n_pass_ewkSsLoose     , w, cw)<<endl;
    cout<<"ewkSsLea         : "<<lcpet(n_pass_ewkSsLea       , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
  }// end for(w)
}
//-----------------------------------------
float SusySelection::computeChargeFlipProb(LeptonVector &leptons, Met &met,
                                           uint systematic, // DG todo
                                           bool update4mom)
{
    return 0.0;
  // cvl_t &ls = leptons;
  // if(ls.size()<2 || !ls[0] || !ls[1] || !m_chargeFlip) return 0.0;
  // Lepton *l0(ls[0]), *l1(ls[1]);
  // int pdg0(susy::pdgIdFromLep(l0)), pdg1(susy::pdgIdFromLep(l1));
  // TLorentzVector smearedLv0(*l0), smearedLv1(*l1);
  // TVector2 smearedMet(met.lv().Px(), met.lv().Py());
  // int sys(NtSys_NOM==systematic ? 0 : 0);
  // //(DGSys_BKGMETHOD_UP==systematic ? +1 : // DG todo : implement syst
  // // (DGSys_BKGMETHOD_DN==systematic ? -1 : 0)));
  // /*
  // cout<<"OS2SS args: "
  //     <<" event   "<<nt.evt()->event
  //     <<" pdg0 "<<pdg0
  //     <<" lv0 px: "<<smearedLv0.Px()<<" py: "<<smearedLv0.Py()<<" pz: "<<smearedLv0.Pz()
  //     <<" pdg1 "<<pdg1
  //     <<" lv1 px: "<<smearedLv1.Px()<<" py: "<<smearedLv1.Py()<<" pz: "<<smearedLv1.Pz()
  //     <<" met px: "<<smearedMet.Px()<<" py: "<<smearedMet.Py()
  //     <<endl;
  // */
  // m_chargeFlip->setSeed(nt.evt()->event);
  // float flipProb(m_chargeFlip->OS2SS(pdg0, &smearedLv0, pdg1, &smearedLv1, &smearedMet, sys));
  // float overlapFrac(m_chargeFlip->overlapFrac().first);
  // if(update4mom) {
  //   m_unsmeared_lv0 = (*l0);
  //   m_unsmeared_lv1 = (*l1);
  //   m_unsmeared_met = met;
  //   l0->SetPtEtaPhiM(smearedLv0.Pt(), smearedLv0.Eta(), smearedLv0.Phi(), smearedLv0.M());
  //   l1->SetPtEtaPhiM(smearedLv1.Pt(), smearedLv1.Eta(), smearedLv1.Phi(), smearedLv1.M());
  //   met.Et = smearedMet.Mod();
  //   met.phi = smearedMet.Phi();
  // }
  // return flipProb*overlapFrac;
}
//-----------------------------------------
susy::wh::Chan SusySelection::getChan(const LeptonVector& leps)
{
  uint ie = 0;
  uint im = 0;
  for(uint i=0; i<leps.size(); ++i){
    if( leps.at(i)->isEle() ) ie++;
    else if( leps.at(i)->isMu() ) im++;
  }
  if( ie == 2 && im == 0 ) return susy::wh::Ch_ee;
  if( ie == 1 && im == 1 ) return susy::wh::Ch_em;
  if( ie == 0 && im == 2 ) return susy::wh::Ch_mm;
  cout<<"Not ee/mm/em... Number Electrons: "<<ie<<" Number Muons: "<<im<<endl;
  return susy::wh::Ch_N; // not in range
}
//-----------------------------------------
SsPassFlags SusySelection::assignNjetFlags(const JetVector& jets, SsPassFlags f)
{
  int njCl = numberOfCLJets(jets);
  int njB  = numberOfCBJets(jets);
  int njF  = numberOfFJets (jets);
  f.bjveto = njB  == 0;
  f.fjveto = njF  == 0;
  f.ge1j   = njCl >= 1;
  f.eq1j   = njCl == 1;
  f.ge2j   = njCl >= 2;
  return f;
}
//-----------------------------------------
bool SusySelection::passThirdLeptonVeto(const Susy::Lepton* l0, const Susy::Lepton* l1, const LeptonVector& otherLeptons, bool verbose)
{
    struct LepPair {
        const Susy::Lepton *signal, *other;
        LepPair(const Susy::Lepton *s, const Susy::Lepton *o) : signal(s), other(o) { assert(s!=0 && o!=0); }
        static const bool unbiasedD0Z0() { return true; }
        static bool lepIsFromPv(const Susy::Lepton* l, float maxD0sig, float maxZ0sinTheta) {
            return (fabs(l->d0Sig(unbiasedD0Z0())) < maxD0sig && fabs(l->z0SinTheta(unbiasedD0Z0())) < maxZ0sinTheta);
        }
        bool otherLepIsElFromPv() { return other->isEle() && lepIsFromPv(other, ELECTRON_D0SIG_CUT_WH, ELECTRON_Z0_SINTHETA_CUT); }
        bool otherLepIsMuFromPv() { return other->isMu()  && lepIsFromPv(other, MUON_D0SIG_CUT,        MUON_Z0_SINTHETA_CUT); }
        bool otherLepIsFromPv() { return otherLepIsElFromPv() || otherLepIsMuFromPv(); }
        bool haveOppositeSign() { return (signal->q * other->q) < 0; }
        bool haveSameFlavor() { return (signal->isMu() && other->isMu()) || (signal->isEle() && other->isEle()); }
        float dR() { return signal->DeltaR(*other); }
        bool areSeparated() { return dR() > 0.05; }
        bool isZcandidate() { return haveOppositeSign() && haveSameFlavor() && areSeparated() && otherLepIsFromPv(); }
        float m() { return (*signal + *other).M(); }
        bool isInZwindow(float maxDelta) { const float mz(91.2); return isZcandidate() && abs(m() - mz) < maxDelta; }
    };
    float maxDeltaMz(20.0);
    size_t nCandInWindow(0);
    for(size_t i=0; i<otherLeptons.size(); ++i) {
        const Susy::Lepton* ol = otherLeptons[i];
        if(verbose) cout<<" lep["<<i<<"] "<<(ol->isEle() ? "E" : ol->isMu() ? "M" : "?")<<" pt "<<ol->Pt()<<endl;
        LepPair ll0(l0, ol), ll1(l1, ol);
        if(ll0.isZcandidate() && ll0.isInZwindow(maxDeltaMz)) {
            nCandInWindow++;
            if(verbose) cout<<" Z candidate: m_ll "<<ll0.m()<<" l0 pt "<<l0->Pt()<<" ol pt "<<ol->Pt()<<" dR "<<ll0.dR()<<endl;
        }
        if(ll1.isZcandidate() && ll1.isInZwindow(maxDeltaMz)) {
            nCandInWindow++;
            if(verbose) cout<<" Z candidate: m_ll "<<ll1.m()<<" l1 pt "<<l1->Pt()<<" ol pt "<<ol->Pt()<<" dR "<<ll1.dR()<<endl;
        }
    }
    if(verbose && nCandInWindow==0) cout<<" No Z candidate"<<endl;
    return nCandInWindow == 0;
}
//-----------------------------------------
void SusySelection::resetAllCounters()
{
  for(int w=0; w<kWeightTypesN; ++w){// Loop over weight types
    n_readin          [w] = 0;
    n_pass_Grl        [w] = 0;
    n_pass_LarErr     [w] = 0;
    n_pass_TileErr    [w] = 0;
    n_pass_TTCVeto    [w] = 0;
    n_pass_GoodVtx    [w] = 0;
    n_pass_TileTrip   [w] = 0;
    n_pass_LAr        [w] = 0;
    n_pass_BadJet     [w] = 0;
    n_pass_FEBCut     [w] = 0;
    n_pass_BadMuon    [w] = 0;
    n_pass_Cosmic     [w] = 0;
    n_pass_hfor       [w] = 0;
    n_pass_HttVeto    [w] = 0;
    n_pass_ge2l       [w] = 0;
    n_pass_eq2l       [w] = 0;
    n_pass_mll        [w] = 0;
    n_pass_signalLep  [w] = 0;
    for(int i=0; i<ET_N; ++i){ // loop over weight x channel.
      // per-SR counters
      n_pass_SR6sign[i][w] = n_pass_SR6flav[i][w] = n_pass_SR6metr[i][w] = 0;
      n_pass_SR6ge1j[i][w] = n_pass_SR6ge2j[i][w] = n_pass_SR6eq2j[i][w] = 0;
      n_pass_SR6eq2jNfv[i][w] = n_pass_SR6ge2jNfv[i][w] = n_pass_SR6[i][w] = 0;
      n_pass_SR6DrllMax     [i][w] = n_pass_SR6PtllMin     [i][w] = 0;
      n_pass_SR6MllMax      [i][w] = n_pass_SR6METRel      [i][w] = 0;
      n_pass_SR6MtLlmetMin  [i][w] = n_pass_SR6MtMinlmetMin[i][w] = 0;
      n_pass_SR6ZtautauVeto [i][w] = 0;

      n_pass_SR7sign[i][w] = n_pass_SR7flav[i][w] = n_pass_SR7metr[i][w] = 0;
      n_pass_SR7ge1j[i][w] = n_pass_SR7ge2j[i][w] = n_pass_SR7eq2j[i][w] = 0;
      n_pass_SR7eq2jNfv[i][w] = n_pass_SR7ge2jNfv[i][w] = n_pass_SR7[i][w] = 0;
      n_pass_SR7DrllMax     [i][w] = n_pass_SR7PtllMin     [i][w] = 0;
      n_pass_SR7MllMax      [i][w] = n_pass_SR7METRel      [i][w] = 0;
      n_pass_SR7MtLlmetMin  [i][w] = n_pass_SR7MtMinlmetMin[i][w] = 0;
      n_pass_SR7ZtautauVeto [i][w] = 0;

      n_pass_SR8sign[i][w] = n_pass_SR8flav[i][w] = n_pass_SR8metr[i][w] = 0;
      n_pass_SR8ge1j[i][w] = n_pass_SR8ge2j[i][w] = n_pass_SR8eq2j[i][w] = 0;
      n_pass_SR8eq2jNfv[i][w] = n_pass_SR8ge2jNfv[i][w] = n_pass_SR8[i][w] = 0;

      n_pass_SR9sign[i][w] = n_pass_SR9flav[i][w] = n_pass_SR9metr[i][w] = 0;
      n_pass_SR9ge1j[i][w] = n_pass_SR9ge2j[i][w] = n_pass_SR9eq2j[i][w] = 0;
      n_pass_SR9eq2jNfv[i][w] = n_pass_SR9ge2jNfv[i][w] = n_pass_SR9[i][w] = 0;

      n_pass_flavor         [i][w] = 0;
      n_pass_os             [i][w] = 0;
      n_pass_ss             [i][w] = 0;
      n_pass_tr2L           [i][w] = 0;
      n_pass_tr2LMatch      [i][w] = 0;
      n_pass_mcTrue2l       [i][w] = 0;
      n_pass_category       [i][w] = 0;
      n_pass_nSigLep        [i][w] = 0;
      n_pass_tauVeto        [i][w] = 0;
      n_pass_mllMin         [i][w] = 0;
      n_pass_muIso          [i][w] = 0;
      n_pass_elD0Sig        [i][w] = 0;
      n_pass_fjVeto         [i][w] = 0;
      n_pass_bjVeto         [i][w] = 0;
      n_pass_ge1j           [i][w] = 0;
      n_pass_eq1j           [i][w] = 0;
      n_pass_ge2j           [i][w] = 0;
      n_pass_lepPt          [i][w] = 0;
      n_pass_mllZveto       [i][w] = 0;
      n_pass_mWwt           [i][w] = 0;
      n_pass_ht             [i][w] = 0;
      n_pass_metRel         [i][w] = 0;
      n_pass_3rdLep         [i][w] = 0;
      n_pass_eq1jlepPt      [i][w] = 0;
      n_pass_eq1jmllZveto   [i][w] = 0;
      n_pass_eq1jDetall     [i][w] = 0;
      n_pass_eq1jMtMax      [i][w] = 0;
      n_pass_eq1jMlj        [i][w] = 0;
      n_pass_eq1jht         [i][w] = 0;
      n_pass_eq1jmetRel     [i][w] = 0;
      n_pass_eq1jmWwt       [i][w] = 0;
      n_pass_ge2jlepPt      [i][w] = 0;
      n_pass_ge2jmllZveto   [i][w] = 0;
      n_pass_ge2jDetall     [i][w] = 0;
      n_pass_ge2jMtMax      [i][w] = 0;
      n_pass_ge2jMljj       [i][w] = 0;
      n_pass_ge2jht         [i][w] = 0;
      n_pass_ge2jmetRel     [i][w] = 0;
      n_pass_ge2jmWwt       [i][w] = 0;
      n_pass_cr1jzv         [i][w] = 0;
      n_pass_cr1jfake       [i][w] = 0;
      n_pass_cr1jzvfake     [i][w] = 0;
      n_pass_cr2jzv         [i][w] = 0;
      n_pass_cr2jfake       [i][w] = 0;
      n_pass_cr2jzvfake     [i][w] = 0;

      n_pass_ewkSs          [i][w] = 0;
      n_pass_ewkSsLoose     [i][w] = 0;
      n_pass_ewkSsLea       [i][w] = 0;
    } // end for(i)
  } // end for(w)
}
//-----------------------------------------
void SusySelection::initChargeFlipTool()
{
  // char* rcdir = getenv("ROOTCOREDIR");
  // if(!rcdir){
  //   if(m_dbg) cout<<"invalid ROOTCOREDIR, cannot initialize chargeFlipTool"<<endl;
  //   return;
  // }
  // string chargeFlipInput(rcdir);
  // chargeFlipInput += "/../ChargeFlip/data/d0_chargeflip_map.root";
  // m_chargeFlip = new chargeFlip(chargeFlipInput);
  // if(m_dbg) m_chargeFlip->printSettings();
}
//-----------------------------------------
LeptonVector SusySelection::getAnyElOrMu(SusyNtObject &susyNt/*, SusyNtSys sys*/)
{
    // DG 2013-12-02:
    // todo1 : re-implement with std algo
    // todo2 : re-implement with syst
    float minPt = 6.0;
    LeptonVector leptons;
    for(uint ie=0; ie<susyNt.ele()->size(); ++ie){
        if(Electron* e = & susyNt.ele()->at(ie)){ //e->setState(sys);
            if(e->Pt()>minPt) leptons.push_back(static_cast<Lepton*>(e));
        }
    }
    for(uint im=0; im<susyNt.muo()->size(); ++im){
        if(Muon* m = & susyNt.muo()->at(im)){ //m->setState(sys);
            if(m->Pt(),minPt) leptons.push_back(static_cast<Lepton*>(m));
        }
    }
    return leptons;
}
//-----------------------------------------
bool SusySelection::passCrWhZVfakeEe(const susy::wh::kin::DilepVars &v)
{
    return (v.isEe
            && fabs(v.mll - 91.2)>10.0
            && v.metrel > 40.0
            && ((v.numCentralLightJets==1 && v.mlj  > 90.0)
                ||
                (v.numCentralLightJets >1 && v.mljj >120.0)));
}
//-----------------------------------------
bool SusySelection::passCrWhZVfakeEm(const susy::wh::kin::DilepVars &v)
{
    return (v.isEm
            && v.pt0 > 30.0
            && v.pt1 > 30.0
            && ((v.numCentralLightJets==1 && v.mlj  > 90.0)
                ||
                (v.numCentralLightJets >1 && v.mljj >120.0)));
}
//-----------------------------------------
bool SusySelection::passCrWhfakeEm  (const susy::wh::kin::DilepVars &v)
{
    return (v.isEm
            && v.pt0 > 30.0
            && v.pt1 < 30.0 // orthogonal to WhZVfake1jem
            && ((v.numCentralLightJets==1 && v.mlj  > 90.0)
                ||
                (v.numCentralLightJets >1 && v.mljj >120.0)));
}
//-----------------------------------------
bool SusySelection::passCrWhZVMm    (const susy::wh::kin::DilepVars &v)
{
    return (v.isMm
            && v.pt0 > 30.0
            && v.pt1 > 30.0
            && ((v.numCentralLightJets==1 && v.mlj  > 90.0)
                ||
                (v.numCentralLightJets >1 && v.mljj >120.0)));
}
//-----------------------------------------
bool SusySelection::passCrWhfakeMm  (const susy::wh::kin::DilepVars &v)
{
    return (v.isMm
            // && v.pt0 > 30.0 // ??
            && v.pt1 < 30.0
            && ((v.numCentralLightJets==1 && v.mlj  > 90.0)
                ||
                (v.numCentralLightJets >1 && v.mljj >120.0)));
}
//-----------------------------------------
bool SusySelection::passCrWhZVfake(const susy::wh::kin::DilepVars &v)
{
    return (SusySelection::passCrWhZVfakeEe(v)
            ||
            SusySelection::passCrWhZVfakeEm(v));
}
//-----------------------------------------
bool SusySelection::passCrWhfake(const susy::wh::kin::DilepVars &v)
{
    return (SusySelection::passCrWhfakeEm(v)
            ||
            SusySelection::passCrWhfakeMm(v));
}
//-----------------------------------------
bool SusySelection::passCrWhZV(const susy::wh::kin::DilepVars &v)
{
    return SusySelection::passCrWhZVMm(v);
}
//-----------------------------------------
DiLepEvtType vars2type(const susy::wh::kin::DilepVars &v)
{
    DiLepEvtType et = ET_Unknown;
    if     (v.isEe) et = ET_ee;
    else if(v.isEm) et = ET_em;
    else if(v.isMm) et = ET_mm;
    return et;
}
//-----------------------------------------
bool SusySelection::passAndIncrementCrWhZVfake(const susy::wh::kin::DilepVars &v)
{
    bool pass(passCrWhZVfake(v));
    if(pass) {
        DiLepEvtType ll = vars2type(v);
        float* counters = v.numCentralLightJets==1 ? n_pass_cr1jzvfake[ll] : n_pass_cr2jzvfake[ll];
        increment(counters, m_weightComponents);
    }
    return pass;
}
//-----------------------------------------
bool SusySelection::passAndIncrementCrWhfake(const susy::wh::kin::DilepVars &v)
{
    bool pass(passCrWhfake(v));
    if(pass) {
        DiLepEvtType ll = vars2type(v);
        float* counters = v.numCentralLightJets==1 ? n_pass_cr1jfake[ll] : n_pass_cr2jfake[ll];
        increment(counters, m_weightComponents);
    }
    return pass;
}
//-----------------------------------------
bool SusySelection::passAndIncrementCrWhZV(const susy::wh::kin::DilepVars &v)
{
    bool pass(passCrWhZV(v));
    if(pass) {
        DiLepEvtType ll = vars2type(v);
        float* counters = v.numCentralLightJets==1 ? n_pass_cr1jzv[ll] : n_pass_cr2jzv[ll];
        increment(counters, m_weightComponents);
    }
    return pass;
}
//-----------------------------------------
bool SusySelection::passSrWh1j(const susy::wh::kin::DilepVars &v, SsPassFlags &f)
{
    bool notApplied(true), pass(false);
    if(v.numCentralLightJets==1){
        if(v.isMm) {
            f.lepPt   = (v.pt0 > 30.0 && v.pt1 > 20.0);
            f.zllVeto = notApplied;
            f.dEtall  = v.detall <   1.5;
            f.maxMt   = v.mtmax()> 100.0;
            f.mljj    = v.mlj    <  90.0;
            f.ht      = v.ht     > 200.0;
            f.metrel  = notApplied;
            f.mtllmet = notApplied;
        } else if(v.isEm) {
            f.lepPt   = (v.pt0 > 30.0 && v.pt1 > 30.0);
            f.zllVeto = notApplied;
            f.dEtall  = v.detall <   1.5;
            f.maxMt   = v.mtmax()> 110.0;
            f.mljj    = v.mlj    <  90.0;
            f.ht      = notApplied;
            f.metrel  = notApplied;
            f.mtllmet = v.mtllmet> 120.0;
        } else if(v.isEe) {
            f.lepPt   = (v.pt0 > 30.0 && v.pt1 > 20.0);
            f.zllVeto = fabs(v.mll-91.2) > 10.0;
            f.dEtall  = notApplied;
            f.maxMt   = notApplied;
            f.mljj    = v.mlj    <  90.0;
            f.ht      = v.ht     > 200.0;
            f.metrel  = v.metrel >  55.0;
            f.mtllmet = notApplied;
        }
        pass = f.lepPt && f.zllVeto && f.dEtall && f.maxMt && f.ht && f.metrel && f.mljj && f.mtllmet;
    }
    return pass;
}
//-----------------------------------------
bool SusySelection::passSrWh1j(const susy::wh::kin::DilepVars &v)
{
    SsPassFlags f;
    return SusySelection::passSrWh1j(v, f);
}
//-----------------------------------------
bool SusySelection::passSrWh2j(const susy::wh::kin::DilepVars &v, SsPassFlags &f)
{
    bool notApplied(true), pass(false);
    if(v.numCentralLightJets>1 && v.numCentralLightJets<4){
        if(v.isMm) {
            f.lepPt = (v.pt0 > 30.0 && v.pt1 > 30.0);
            f.zllVeto = notApplied;
            f.dEtall  = v.detall <   1.5;
            f.maxMt   = notApplied;
            f.mljj    = v.mljj   < 120.0;
            f.ht      = v.ht     > 220.0;
            f.metrel  = notApplied;
            f.mtllmet = notApplied;
        } else if(v.isEm) {
            f.lepPt   = (v.pt0 > 30.0 && v.pt1 > 30.0);
            f.zllVeto = notApplied;
            f.dEtall  = v.detall <   1.5;
            f.maxMt   = notApplied;
            f.mljj    = v.mljj   < 120.0;
            f.ht      = notApplied;
            f.metrel  = notApplied;
            f.mtllmet = v.mtllmet> 110.0;
        } else if(v.isEe) {
            f.lepPt   = (v.pt0 > 30.0 && v.pt1 > 20.0);
            f.zllVeto = fabs(v.mll-91.2) > 10.0;
            f.dEtall  = notApplied;
            f.maxMt   = v.mtmax()> 100.0;
            f.mljj    = v.mljj   < 120.0;
            f.ht      = notApplied;
            f.metrel  = v.metrel >  30.0;
            f.mtllmet = notApplied;
        }
    pass = f.lepPt && f.zllVeto && f.dEtall && f.maxMt && f.mljj && f.ht && f.metrel && f.mtllmet;
    }
    return pass;
}
//-----------------------------------------
bool SusySelection::passSrWh2j(const susy::wh::kin::DilepVars &v)
{
    SsPassFlags f;
    return SusySelection::passSrWh2j(v, f);
}
//-----------------------------------------
bool SusySelection::passSrWhNoMlj(const susy::wh::kin::DilepVars &v)
{
    SsPassFlags f1j, f2j;
    passSrWh1j(v, f1j);
    passSrWh2j(v, f2j);
    f1j.mljj = true;
    f2j.mljj = true;
    bool pass1j(f1j.lepPt && f1j.zllVeto && f1j.dEtall && f1j.maxMt && f1j.mljj && f1j.ht && f1j.metrel && f1j.mtllmet);
    bool pass2j(f2j.lepPt && f2j.zllVeto && f2j.dEtall && f2j.maxMt && f2j.mljj && f2j.ht && f2j.metrel && f2j.mtllmet);
    return (pass1j || pass2j);
}
//-----------------------------------------
bool SusySelection::passSrRazor0jet(const LeptonVector &leptons, const JetVector& jets, const Met &met)
{
    bool pass=false;
    if(leptons.size()<2) return pass;
    size_t num_central_light_jets = numberOfCLJets(jets);
    size_t num_central_cb_jets = numberOfCBJets(jets);
    size_t num_forward_jets = numberOfFJets(jets);
    Susy::Lepton &l0 = *leptons[0];
    Susy::Lepton &l1 = *leptons[1];
    bool ee(l0.isEle() && l1.isEle());
    bool mumu(l0.isMu() && l1.isMu());
    bool emu((l0.isEle() && l1.isMu()) || (l0.isMu() && l1.isEle()));
    double mll((l0+l1).M());
    double dphi_ll_vBetaT(0.0), mDeltaR(0.0);
    SusySelection::computeRazor(leptons, jets, met, dphi_ll_vBetaT, mDeltaR);

    if(ee){
        pass =(abs(mll-91.2) > 10.0 &&
               // num_forward_jets==0 &&
               // num_central_cb_jets==0 &&
               // num_central_light_jets==0 &&
               l0.Pt() > 20.0 &&
               l1.Pt() > 20.0 //&&
               // mDeltaR > 20.0
               // mDeltaR>150.0
            );
    } else if(mumu) {
        pass =(abs(mll-91.2) > 10.0 &&
               // num_forward_jets==0 &&
               // num_central_cb_jets==0 &&
               // num_central_light_jets==0 &&
               l0.Pt() > 20.0 &&
               l1.Pt() > 20.0 //&&
               // mDeltaR > 20.0
               // mDeltaR > 150.0
            );
    } else if(emu) {
        pass =(// num_forward_jets==0 &&
               // num_central_cb_jets==0 &&
               // num_central_light_jets==0 &&
               l0.Pt() > 20.0 &&
               l1.Pt() > 20.0 //&&
               // mDeltaR > 20.0
               // mDeltaR > 150.0
            );
    }
//    cout<<"SusySelection::passSrRazor0jet "<<(ee ? "ee" : mumu ? "mumu" : "emu")<<" "<<pass<<endl;
    return pass;
}
//-----------------------------------------
bool SusySelection::passSrRazor1jet(const LeptonVector &leptons, const JetVector& jets, const Met &met)
{
    bool pass=false;
    if(leptons.size()<2) return pass;
    size_t num_central_light_jets = numberOfCLJets(jets);
    size_t num_central_cb_jets = numberOfCBJets(jets);
    size_t num_forward_jets = numberOfFJets(jets);
    Susy::Lepton &l0 = *leptons[0];
    Susy::Lepton &l1 = *leptons[1];
    bool ee(l0.isEle() && l1.isEle());
    bool mumu(l0.isMu() && l1.isMu());
    bool emu((l0.isEle() && l1.isMu()) || (l0.isMu() && l1.isEle()));
    double mll((l0+l1).M());
    double dphi_ll_vBetaT(0.0), mDeltaR(0.0);
    SusySelection::computeRazor(leptons, jets, met, dphi_ll_vBetaT, mDeltaR);
    double r2 = (met.Et / (met.Et + l0.Pt() + l1.Pt()));
    TLorentzVector j0 = m_signalJets2Lep.size()>0 ? *m_signalJets2Lep[0] : TLorentzVector();
    if(ee){
        pass = (num_central_light_jets==1 &&
                num_central_cb_jets==0 &&
                num_forward_jets==0 &&
                abs(mll-91.2) > 10.0 &&
                mDeltaR > 90.0 &&
                j0.Pt() > 80.0 &&
                dphi_ll_vBetaT > 2.0 &&
                r2 > 0.5);
    } else if(mumu) {
        pass = (num_central_light_jets==1 &&
                num_central_cb_jets==0 &&
                num_forward_jets==0 &&
                abs(mll-91.2) > 10.0 &&
                mDeltaR > 90.0 &&
                j0.Pt() > 80.0 &&
                dphi_ll_vBetaT > 2.0 &&
                r2 > 0.5);
    } else if(emu) {
        pass = (num_central_light_jets==1 &&
                num_central_cb_jets==0 &&
                num_forward_jets==0 &&
                mDeltaR > 90.0 &&
                j0.Pt() > 80.0 &&
                dphi_ll_vBetaT > 2.5 &&
                r2 > 0.7);
    }
//    cout<<"SusySelection::passSrRazor0jet "<<(ee ? "ee" : mumu ? "mumu" : "emu")<<" "<<pass<<endl;
    return pass;
}
//-----------------------------------------

void SusySelection::computeRazor(const LeptonVector &leptons, const JetVector& jets, const Met &met,
                                 double &dphi_ll_vBetaT, double &mDeltaR)
{
    TVector3 vBETA_z, pT_CM, vBETA_T_CMtoR, vBETA_R;
    double SHATR, dphi_LL_vBETA_T, dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T, MDELTAR, costhetaRp1;
    SusyNtTools::superRazor(leptons, &met,
                            vBETA_z,  pT_CM, vBETA_T_CMtoR,  vBETA_R,
                            SHATR, dphi_LL_vBETA_T, dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T, MDELTAR, costhetaRp1);
    dphi_ll_vBetaT = dphi_LL_vBETA_T;
    mDeltaR = MDELTAR;
}
