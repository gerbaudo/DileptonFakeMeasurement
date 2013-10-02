#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <sstream>
#include "TCanvas.h"
#include "SusyTest0/SusySelection.h"
#include "SusyTest0/SusyPlotter.h"

#include "LeptonTruthTools/RecoTruthMatch.h" // provides RecoTruthMatch::
#include "ChargeFlip/chargeFlip.h"
#include "SusyTest0/criteria.h"

using namespace std;
using namespace Susy;

SusySelection::SusySelection() :
  m_susyObj(NULL),
  m_xsReader(NULL),
  m_trigObj(NULL),
  m_useMCTrig(false),
  m_w(1.0),
  m_useXsReader(false),
  m_xsFromReader(-1.0),
  m_qflipProb(0.0)
{
  resetAllCounters();
  setAnaType(Ana_2LepWH);
  setSelectTaus(true);
  initChargeFlipTool();
}
void SusySelection::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  if(m_dbg) cout << "SusySelection::Begin" << endl;
  string period = "Moriond";
  bool useReweightUtils = false;
  m_trigObj = new DilTrigLogic(period, useReweightUtils);
  if(m_useMCTrig) m_trigObj->useMCTrigger();
  if( m_useXsReader ){
    m_xsReader = new XSReader();
    m_xsReader->setDebug(m_dbg);
    m_xsReader->LoadXSInfo();
  } // end if(m_useXsReader)
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
  passSrSs(WH_SRSS1, m_signalLeptons, m_signalTaus, m_signalJets2Lep, m_met, allowQflip);

  return kTRUE;
}
//-----------------------------------------
void SusySelection::Terminate()
{
  SusyNtAna::Terminate();
  if(m_dbg) cout << "SusySelection::Terminate" << endl;
  dumpEventCounters();
  if(m_xsReader) delete m_xsReader;
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
bool SusySelection::passSrSs(const WH_SR signalRegion,
                             LeptonVector& leptons,
                             const TauVector& taus,
                             const JetVector& jets,
                             const Met *met,
                             bool allowQflip)
{
  if(leptons.size()<2) return false;
  DiLepEvtType ll(getDiLepEvtType(leptons));
  const DiLepEvtType ee(ET_ee), em(ET_em), me(ET_me), mm(ET_mm);
  WH_SR sr = signalRegion;
  float ptL0Min  = 30;
  float ptL1Min  = (ll==mm ? 0.0 : 20.0);
  float htMin    = 200;
  bool applyMllZveto(ll==ee);
  float mZ0(91.2);
  float loMllZ(applyMllZveto ? mZ0-10. : FLT_MAX);
  float hiMllZ(applyMllZveto ? mZ0+10. : FLT_MIN);
  float mtwwMin = (ll==ee ? 150 : (ll==em || ll==me ? 140 :
                                   (ll==mm ? (sr==WH_SRSS1 ? 100 :
                                              (sr==WH_SRSS2 ? 150 :
                                               (sr==WH_SRSS3 ? 200 :
                                                FLT_MIN))) : FLT_MIN)));
  float metRelMin = (ll==ee ? 50 : (ll==em || ll==me ? 50 :
                                    (ll==mm ? (sr==WH_SRSS4 ? 50 :
                                               FLT_MIN) : FLT_MIN)));
  bool update4mom(true); // charge flip
  bool u4m=update4mom;
  // if(!passEventCleaning()){ return false; }

  // Apply event selection cuts
  bool mc(nt.evt()->isMC), data(!mc);
  const LeptonVector &ls = leptons;
  LeptonVector &ncls = leptons; // non-const leptons: can be modified by qflip
  Met ncmet(*m_met); // non-const met \todo: should modify a non-const input
  const JetVector    &js = jets;
  WeightComponents &wc = m_weightComponents;

  if(true)                                      increment(n_pass_category [ll], wc); else return false;
  if(susy::passNlepMin(ls, 2))                  increment(n_pass_nSigLep  [ll], wc); else return false;
  if(m_signalTaus.size()==0)                    increment(n_pass_tauVeto  [ll], wc); else return false;
  if(passTrig2L     (ls))                       increment(n_pass_tr2L     [ll], wc); else return false;
  if(passTrig2LMatch(ls))                       increment(n_pass_tr2LMatch[ll], wc); else return false;
  if(data || susy::isTrueDilepton(ls))          increment(n_pass_mcTrue2l [ll], wc); else return false;
  bool sameSign = allowQflip ? sameSignOrQflip(ncls, ncmet, ll, u4m, mc) : susy::sameSign(ncls);
  if(sameSign)                                  increment(n_pass_ss       [ll], wc); else return false;
  met = &ncmet; // after qflip, use potentially smeared lep and met
                                                increment(n_pass_muIso    [ll], wc);
                                                increment(n_pass_elD0Sig  [ll], wc);
  if(passfJetVeto  (js))                        increment(n_pass_fjVeto   [ll], wc); else return false;
  if(passbJetVeto  (js))                        increment(n_pass_bjVeto   [ll], wc); else return false;
  if(passge1Jet    (js))                        increment(n_pass_ge1j     [ll], wc); else return false;
  if(susy::pass2LepPt(ncls, ptL0Min, ptL1Min))  increment(n_pass_lepPt    [ll], wc); else return false;
  if(susy::passZllVeto(ncls, loMllZ, hiMllZ))   increment(n_pass_mllZveto [ll], wc); else return false;
  if(susy::passMtLlMetMin(ncls, met, mtwwMin))  increment(n_pass_mWwt     [ll], wc); else return false;
  if(susy::passHtMin(ncls, js, met, htMin))     increment(n_pass_ht       [ll], wc); else return false;
  if(passMetRelMin (met,ncls,js,metRelMin))     increment(n_pass_metRel   [ll], wc); else return false;

  /*
  cout<<""
      <<" ("<<(ll==ee ? "EE" : (ll==mm ? "MM" : "EM"))<<")"
      <<" event   "<<nt.evt()->event
      <<" gen     "<<wc.gen
      <<" pileup  "<<wc.pileup
      <<" norm    "<<wc.norm
      <<" lepSf   "<<wc.lepSf
      <<" btag    "<<wc.btag
      <<" trigger "<<wc.trigger
      <<" qflip   "<<wc.qflip
      <<" all     "<<wc.product()
      <<endl;
  */
  return true;
}
//-----------------------------------------
bool SusySelection::passHfor()
{
  if(nt.evt()->hfor == 4 ) return false;
  return true;
}
//-----------------------------------------
bool SusySelection::passTrig2L(const LeptonVector& leptons)
{
  if(leptons.size() != 2 || !m_trigObj) return false;
  return m_trigObj->passDilEvtTrig(leptons, m_met->Et, nt.evt());
}
//-----------------------------------------
bool SusySelection::passTrig2LMatch(const LeptonVector& leptons)
{
  if(leptons.size() != 2 || !m_trigObj) return false;
  return m_trigObj->passDilTrigMatch(leptons, m_met->Et, nt.evt());
}
//-----------------------------------------
bool SusySelection::passTrig2LwithMatch(const LeptonVector& leptons)
{
  return (passTrig2L(leptons) && passTrig2LwithMatch(leptons));
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
  // if m_useMCTrig, then we are dropping evts with DilTrigLogic::passDil*, not weighting them
  // DG Jun2013 : verify this with Matt & Josephine
  if(!m_useMCTrig){
    if(leptons.size()==2) trigW = m_trigObj->getTriggerWeight(leptons,
                                                              nt.evt()->isMC,
                                                              m_met->Et,
                                                              m_signalJets2Lep.size(),
                                                              nt.evt()->nVtx,
                                                              NtSys_NOM);
    bool twIsInvalid(isnan(trigW) || trigW<0.0);
    assert(!twIsInvalid);
    if(twIsInvalid){
      if(m_dbg)
        cout<<"SusySelection::getTriggerWeight: invalid weigth "<<trigW<<", using 0.0"<<endl;
      trigW = (twIsInvalid ? 0.0 : trigW);
    }
  }
  return trigW;
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
    cout<<"fjVeto           : "<<lcpet(n_pass_fjVeto         , w, cw)<<endl;
    cout<<"bjVeto           : "<<lcpet(n_pass_bjVeto         , w, cw)<<endl;
    cout<<"ge1j             : "<<lcpet(n_pass_ge1j           , w, cw)<<endl;
    cout<<"lepPt            : "<<lcpet(n_pass_lepPt          , w, cw)<<endl;
    cout<<"mllZveto         : "<<lcpet(n_pass_mllZveto       , w, cw)<<endl;
    cout<<"mWwt             : "<<lcpet(n_pass_mWwt           , w, cw)<<endl;
    cout<<"ht               : "<<lcpet(n_pass_ht             , w, cw)<<endl;
    cout<<"metRel           : "<<lcpet(n_pass_metRel         , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
    cout<<"mllMin           : "<<lcpet(n_pass_mllMin         , w, cw)<<endl;
    cout<<"flavor:          : "<<lcpet(n_pass_flavor         , w, cw)<<endl;
    cout<<"OS:              : "<<lcpet(n_pass_os             , w, cw)<<endl;
    cout<<midRule                                                    <<endl;
  }// end for(w)
}
//-----------------------------------------
float SusySelection::getXsFromReader()
{
  if(!m_useXsReader || !m_xsReader) return -1.0;
  bool xsIsNotCached(m_xsFromReader < 0.0); // was initialized to -1
  if(xsIsNotCached){
    int dsid(static_cast<int>(nt.evt()->mcChannel));
    m_xsFromReader = m_xsReader->GetXS(dsid);
    if(m_dbg) cout<<"SusySelection::getXsFromReader: got "<<m_xsFromReader<<" for "<<dsid<<endl;
  }
  return m_xsFromReader;
}
//-----------------------------------------
float SusySelection::computeEventWeightXsFromReader(float lumi)
{
  float defaultXsec = nt.evt()->xsec;
  assert(defaultXsec != 0.0);
  return (getEventWeight(lumi) * getXsFromReader() / defaultXsec);
}
//-----------------------------------------
float SusySelection::computeChargeFlipProb(LeptonVector &leptons, Met &met,
                                           uint systematic, // DG todo
                                           bool update4mom)
{
  cvl_t &ls = leptons;
  if(ls.size()<2 || !ls[0] || !ls[1] || !m_chargeFlip) return 0.0;
  Lepton *l0(ls[0]), *l1(ls[1]);
  int pdg0(susy::pdgIdFromLep(l0)), pdg1(susy::pdgIdFromLep(l1));
  TLorentzVector smearedLv0(*l0), smearedLv1(*l1);
  TVector2 smearedMet(met.lv().Px(), met.lv().Py());
  int sys(NtSys_NOM==systematic ? 0 : 0);
  //(DGSys_BKGMETHOD_UP==systematic ? +1 : // DG todo : implement syst
  // (DGSys_BKGMETHOD_DN==systematic ? -1 : 0)));
  /*
  cout<<"OS2SS args: "
      <<" event   "<<nt.evt()->event
      <<" pdg0 "<<pdg0
      <<" lv0 px: "<<smearedLv0.Px()<<" py: "<<smearedLv0.Py()<<" pz: "<<smearedLv0.Pz()
      <<" pdg1 "<<pdg1
      <<" lv1 px: "<<smearedLv1.Px()<<" py: "<<smearedLv1.Py()<<" pz: "<<smearedLv1.Pz()
      <<" met px: "<<smearedMet.Px()<<" py: "<<smearedMet.Py()
      <<endl;
  */
  m_chargeFlip->setSeed(nt.evt()->event);
  float flipProb(m_chargeFlip->OS2SS(pdg0, &smearedLv0, pdg1, &smearedLv1, &smearedMet, sys));
  float overlapFrac(m_chargeFlip->overlapFrac().first);
  if(update4mom) {
    m_unsmeared_lv0 = (*l0);
    m_unsmeared_lv1 = (*l1);
    m_unsmeared_met = met;
    l0->SetPtEtaPhiM(smearedLv0.Pt(), smearedLv0.Eta(), smearedLv0.Phi(), smearedLv0.M());
    l1->SetPtEtaPhiM(smearedLv1.Pt(), smearedLv1.Eta(), smearedLv1.Phi(), smearedLv1.M());
    met.Et = smearedMet.Mod();
    met.phi = smearedMet.Phi();
  }
  return flipProb*overlapFrac;
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

      n_pass_flavor     [i][w] = 0;
      n_pass_os         [i][w] = 0;
      n_pass_ss         [i][w] = 0;
      n_pass_tr2L       [i][w] = 0;
      n_pass_tr2LMatch  [i][w] = 0;
      n_pass_mcTrue2l   [i][w] = 0;
      n_pass_category   [i][w] = 0;
      n_pass_nSigLep    [i][w] = 0;
      n_pass_tauVeto    [i][w] = 0;
      n_pass_mllMin     [i][w] = 0;
      n_pass_muIso      [i][w] = 0;
      n_pass_elD0Sig    [i][w] = 0;
      n_pass_fjVeto     [i][w] = 0;
      n_pass_bjVeto     [i][w] = 0;
      n_pass_ge1j       [i][w] = 0;
      n_pass_lepPt      [i][w] = 0;
      n_pass_mllZveto   [i][w] = 0;
      n_pass_mWwt       [i][w] = 0;
      n_pass_ht         [i][w] = 0;
      n_pass_metRel     [i][w] = 0;
    } // end for(i)
  } // end for(w)
}
//-----------------------------------------
void SusySelection::initChargeFlipTool()
{
  char* rcdir = getenv("ROOTCOREDIR");
  if(!rcdir){
    if(m_dbg) cout<<"invalid ROOTCOREDIR, cannot initialize chargeFlipTool"<<endl;
    return;
  }
  string chargeFlipInput(rcdir);
  chargeFlipInput += "/../ChargeFlip/data/d0_chargeflip_map.root";
  m_chargeFlip = new chargeFlip(chargeFlipInput);
  if(m_dbg) m_chargeFlip->printSettings();
}
//-----------------------------------------
