#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <sstream>
#include "TCanvas.h"
#include "SusyTest0/SusySelection.h"
#include "SusyTest0/SusyPlotter.h"

#include "Mt2/mt2_bisect.h"
#include "LeptonTruthTools/RecoTruthMatch.h" // provides RecoTruthMatch::
#include "SusyNtuple/WhTruthExtractor.h"
#include "ChargeFlip/chargeFlip.h"

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
  m_ET(ET_Unknown),
  m_qflipProb(0.0)
{
  resetAllCounters();
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
  GetEntry(entry);
  clearObjects();
  m_ET = ET_Unknown;
  m_chainEntry++;
  increment(n_readin);
  if(m_dbg || m_chainEntry%50000==0)
  {
    cout<<"****"
        <<" Processing entry "<<setw(6)<<m_chainEntry
        <<" run "             <<setw(6)<<nt.evt()->run
        <<" event "           <<setw(7)<<nt.evt()->event
        <<" ****"<<endl;
  }
  bool removeLepsFromIso(false);
  selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium);
  if(!selectEvent()) return kTRUE;

  m_ET = getDiLepEvtType(m_baseLeptons);
  passSrSs(m_ET, WH_SRSS1, m_signalLeptons, m_signalTaus, m_signalJets2Lep, m_met);

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
void SusySelection::increment(float counters[],
                              bool includeLepSF, bool includeBtag)
{
  float eventWeight = nt.evt()->w;
  float pileup = nt.evt()->wPileup;
  float lepSf = (includeLepSF ?
                 m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF
                 : 1.0);
  float btag = includeBtag ? getBTagWeight(nt.evt()) : 1.0;
  bool includeTrig(m_baseLeptons.size() == 2 && nt.evt()->isMC); //DG should be arg?
  float trig = (includeTrig ?
                m_trigObj->getTriggerWeight(m_baseLeptons,
                                            nt.evt()->isMC,
                                            m_met->Et,
                                            m_signalJets2Lep.size(),
                                            nt.evt()->nVtx,
                                            NtSys_NOM)
                : 1.0);
  float all = getEventWeightFixed(nt.evt()->mcChannel,LUMI_A_L) *btag*trig*lepSf;
  counters[WT_Raw ] += 1.0;
  counters[WT_Evt ] += eventWeight;
  counters[WT_PU  ] += eventWeight * pileup;
  counters[WT_LSF ] += eventWeight * lepSf;
  counters[WT_Btag] += eventWeight * btag;
  counters[WT_Trig] += eventWeight * trig;
  counters[WT_All ] += all;
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
  if(passGRL        (flag           ))  { increment(n_pass_Grl     );} else { return false; }
  if(passLarErr     (flag           ))  { increment(n_pass_LarErr  );} else { return false; }
  if(passTileErr    (flag           ))  { increment(n_pass_TileErr );} else { return false; }
  if(passTTCVeto    (flag           ))  { increment(n_pass_TTCVeto );} else { return false; }
  if(passGoodVtx    (flag           ))  { increment(n_pass_GoodVtx );} else { return false; }
  if(passTileTripCut(flag           ))  { increment(n_pass_TileTrip);} else { return false; }
  if(passLAr        (flag           ))  { increment(n_pass_LAr     );} else { return false; }
  if(!hasBadJet     (jets           ))  { increment(n_pass_BadJet  );} else { return false; }
  if(passDeadRegions(pjets,met,run,mc)) { increment(n_pass_FEBCut  );} else { return false; }
  if(!hasBadMuon    (m_preMuons     ))  { increment(n_pass_BadMuon );} else { return false; }
  if(!hasCosmicMuon (m_baseMuons    ))  { increment(n_pass_Cosmic  );} else { return false; }
  if(passHfor       (               ))  { increment(n_pass_hfor    );} else { return false; }
  //if(passHtautauVeto(hdec)) { increment(n_pass_HttVeto ); } else { return false; }
  if(bleps.size() >= 2               )  { increment(n_pass_ge2l    );} else { return false; }
  if(bleps.size() == 2               )  { increment(n_pass_eq2l    );} else { return false; }
  if(passMllMin(bleps, mllMin       ))  { increment(n_pass_mll     );} else { return false; }
  return true;
}
//-----------------------------------------
bool SusySelection::passSR6base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return oppositeSign(leptons) && sameFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR7base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return oppositeSign(leptons) && oppositeFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR8base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return sameSign(leptons) && sameFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR9base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return sameSign(leptons) && oppositeFlavor(leptons);
}
//-----------------------------------------
bool SusySelection::passSR6(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  bool lepSf(true), bSf(true);
  if (oppositeSign(leptons)) { if (count) increment(n_pass_SR6sign[m_ET],lepSf, bSf); }
  else return false;
  if (sameFlavor(leptons))   { if (count) increment(n_pass_SR6flav[m_ET],lepSf, bSf); }
  else return false;
  if (count && passge1Jet(jets) )         increment(n_pass_SR6ge1j[m_ET],lepSf, bSf);
  if (passge2Jet(jets))      { if (count) increment(n_pass_SR6ge2j[m_ET],lepSf, bSf); }
  else return false;
  if( passDrllMax(leptons) )               {if(count) increment(n_pass_SR6DrllMax     [m_ET], lepSf, bSf);}
  else return false;
  if( passPtllMin(leptons) )               {if(count) increment(n_pass_SR6PtllMin     [m_ET], lepSf, bSf);}
  else return false;
  if( passMllMax(leptons) )                {if(count) increment(n_pass_SR6MllMax      [m_ET], lepSf, bSf);}
  else return false;
  if(passMetRelMin(met,leptons,jets,50.0)) {if(count) increment(n_pass_SR6METRel      [m_ET], lepSf, bSf);}
  else return false;
  if( passMtLlMetMin(leptons, met) )       {if(count) increment(n_pass_SR6MtLlmetMin  [m_ET], lepSf, bSf);}
  else return false;
  if( passMtMinlmetMin(leptons, met) )     {if(count) increment(n_pass_SR6MtMinlmetMin[m_ET], lepSf, bSf);}
  else return false;
  if( passZtautauVeto(leptons, jets, met)) {if(count) increment(n_pass_SR6ZtautauVeto [m_ET], lepSf, bSf);}
  else return false;
  if( count )                               increment(n_pass_SR6[m_ET],lepSf, bSf);
  return true;
}
//-----------------------------------------
bool SusySelection::passSR7(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  bool lepSf(true), bSf(true);
  if( oppositeSign(leptons) ) { if (count) increment(n_pass_SR7sign[m_ET],lepSf, bSf);}
  else return false;
  if( oppositeFlavor(leptons) ) { if(count) increment(n_pass_SR7flav[m_ET],lepSf, bSf);}
  else return false;
  if( count && passge1Jet(jets) )           increment(n_pass_SR7ge1j[m_ET],lepSf, bSf);
  if( passge2Jet(jets)) { if(count)         increment(n_pass_SR7ge2j[m_ET],lepSf, bSf);}
  else return false;
  if( passDrllMax(leptons) )               {if(count) increment(n_pass_SR7DrllMax     [m_ET], lepSf, bSf);}
  else return false;
  if( passPtllMin(leptons) )               {if(count) increment(n_pass_SR7PtllMin     [m_ET], lepSf, bSf);}
  else return false;
  if( passMllMax(leptons) )                {if(count) increment(n_pass_SR7MllMax      [m_ET], lepSf, bSf);}
  else return false;
  if(passMetRelMin(met,leptons,jets,50.0)) {if(count) increment(n_pass_SR7METRel      [m_ET], lepSf, bSf);}
  else return false;
  if( passMtLlMetMin(leptons, met) )       {if(count) increment(n_pass_SR7MtLlmetMin  [m_ET], lepSf, bSf);}
  else return false;
  if( passMtMinlmetMin(leptons, met) )     {if(count) increment(n_pass_SR7MtMinlmetMin[m_ET], lepSf, bSf);}
  else return false;
  if( passZtautauVeto(leptons, jets, met)) {if(count) increment(n_pass_SR7ZtautauVeto [m_ET], lepSf, bSf);}
  else return false;
  //if( !passPtTot      (leptons, jets, met)) return false;
  if( count )                               increment(n_pass_SR7[m_ET],lepSf, bSf);
  return true;
}
//-----------------------------------------
bool SusySelection::passSR8(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  bool lepSf(true), bSf(true);
  if( sameSign(leptons) )                {if( count ) increment(n_pass_SR8sign[m_ET],lepSf, bSf); }
  else return false;
  if( sameFlavor(leptons) )              {if( count ) increment(n_pass_SR8flav[m_ET],lepSf, bSf);}
  else return false;
  if( count && passge1Jet(jets) )                      increment(n_pass_SR8ge1j[m_ET],lepSf, bSf);
  if( passge2Jet(jets))                  {if( count ) increment(n_pass_SR8ge2j[m_ET],lepSf, bSf); }
  else return false;
  if(passMetRelMin(met,leptons,jets,50.)){if( count ) increment(n_pass_SR8metr[m_ET],lepSf, bSf);}
  else return false;
  return true;
}
//-----------------------------------------
bool SusySelection::passSR9(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  bool lepSf(true), bSf(true);
  if( sameSign(leptons) )                {if( count ) increment(n_pass_SR9sign[m_ET],lepSf, bSf);}
  else return false;
  if( oppositeFlavor(leptons) )          {if( count ) increment(n_pass_SR9flav[m_ET],lepSf, bSf);}
  else return false;
  if( count && passge1Jet(jets) )                     increment(n_pass_SR9ge1j[m_ET],lepSf, bSf);
  if( passge2Jet(jets))                  {if( count ) increment(n_pass_SR9ge2j[m_ET],lepSf, bSf);}
  else return false;
  if(passMetRelMin(met,leptons,jets,50.)){if( count ) increment(n_pass_SR9metr[m_ET],lepSf, bSf);}
  else return false;
  return true;
}
//-----------------------------------------
bool SusySelection::passSrSsBase()
{
  return false;
}
//-----------------------------------------
bool SusySelection::passSrSs(const DiLepEvtType eventType,
                             const WH_SR signalRegion,
                             LeptonVector& leptons,
                             const TauVector& taus,
                             const JetVector& jets,
                             const Met *met)
{
  DiLepEvtType ll = eventType;
  const DiLepEvtType ee(ET_ee), em(ET_em), me(ET_me), mm(ET_mm);
  WH_SR sr = signalRegion;
  float muIsoMax = 0.1;
  float ptL0Min  = 30;
  float ptL1Min  = (ll==mm ? 0.0 : 20.0);
  float htMin    = 200;
  float d0SMax   = (ll==ee || ll==em || ll==me ?   3 : FLT_MAX);
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
  bool lepSf(true), bSf(true);
  bool update4mom(true); // charge flip
  bool u4m=update4mom;
  // if(!passEventCleaning()){ return false; }

  // Apply event selection cuts
  bool mc(nt.evt()->isMC), data(!mc);
  const LeptonVector &ls = leptons;
  LeptonVector &ncls = leptons; // non-const leptons: can be modified by qflip
  Met ncmet(*m_met); // non-const met
  const JetVector    &js = jets;
  if(true)                                      increment(n_pass_category [ll], lepSf, bSf); else return false;
  if(passNlepMin   (ls, 2))                     increment(n_pass_nSigLep  [ll], lepSf, bSf); else return false;
  if(m_signalTaus.size()==0)                    increment(n_pass_tauVeto  [ll], lepSf, bSf); else return false;
  if(passTrig2L     (ls))                       increment(n_pass_tr2L     [ll], lepSf, bSf); else return false;
  if(passTrig2LMatch(ls))                       increment(n_pass_tr2LMatch[ll], lepSf, bSf); else return false;
  if(data || isTrueDilepton(ls))                increment(n_pass_mcTrue2l [ll], lepSf, bSf); else return false;
  if(sameSignOrQflip(ncls, ncmet, ll, u4m, mc)) increment(n_pass_ss       [ll], lepSf, bSf); else return false;
  met = &ncmet; // after qflip, use potentially smeared lep and met
  if(passMuonRelIso(ncls, muIsoMax))            increment(n_pass_muIso    [ll], lepSf, bSf); else return false;
  if(passEleD0S    (ncls, d0SMax))              increment(n_pass_elD0Sig  [ll], lepSf, bSf); else return false;
  if(passfJetVeto  (js))                        increment(n_pass_fjVeto   [ll], lepSf, bSf); else return false;
  if(passbJetVeto  (js))                        increment(n_pass_bjVeto   [ll], lepSf, bSf); else return false;
  if(passge1Jet    (js))                        increment(n_pass_ge1j     [ll], lepSf, bSf); else return false;
  if(pass2LepPt    (ncls, ptL0Min, ptL1Min))    increment(n_pass_lepPt    [ll], lepSf, bSf); else return false;
  if(passZllVeto   (ncls, loMllZ, hiMllZ))      increment(n_pass_mllZveto [ll], lepSf, bSf); else return false;
  if(passMtLlMetMin(ncls, met, mtwwMin))        increment(n_pass_mWwt     [ll], lepSf, bSf); else return false;
  if(passHtMin     (ncls, js, met, htMin))      increment(n_pass_ht       [ll], lepSf, bSf); else return false;
  if(passMetRelMin (met,ncls,js,metRelMin))     increment(n_pass_metRel   [ll], lepSf, bSf); else return false;
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
bool SusySelection::sameFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return (leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
//-----------------------------------------
bool SusySelection::oppositeFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return !(leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
//-----------------------------------------
bool SusySelection::sameSign(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return leptons.at(0)->q * leptons.at(1)->q > 0;
}
//-----------------------------------------
bool SusySelection::sameSignOrQflip(LeptonVector& leptons, Met &met,
                                    const DiLepEvtType eventType,
                                    bool update4mom, bool isMc)
{
  bool isSS(sameSign(leptons));
  if(!isMc) return isSS;
  bool isOS(!isSS);
  bool canBeQflip(isOS && (eventType==ET_ee || eventType==ET_em || eventType==ET_me));
  if (!isSS && !canBeQflip) return false;
  if(canBeQflip){
    uint systematic=NtSys_NOM; // DG sys todo
    m_qflipProb = computeChargeFlipProb(leptons, met, systematic, update4mom);
  }
  return true;
}
//-----------------------------------------
bool SusySelection::oppositeSign(const LeptonVector& leptons)
{
  return !(sameSign(leptons));
}
//-----------------------------------------
bool SusySelection::passHtautauVeto(int hdecay)
{
  return (hdecay!=WhTruthExtractor::kPtauAtau);
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
//-----------------------------------------
bool SusySelection::passdPhi(TLorentzVector v0, TLorentzVector v1, float cut)
{
  return v0.DeltaPhi(v1) > cut;
}
//-----------------------------------------
bool SusySelection::passMtLlMetMin(const LeptonVector& l, const Met* met,
                                   float minVal)
{
  if(l.size() < 2 || !l[0] || !l[1]) return false;
  TLorentzVector ll = (*l[0] + *l[1]);
  float mww=SusyPlotter::transverseMass(ll, met->lv());
  return (minVal < mww);
}
//----------------------------------------------------------
bool SusySelection::passMtMinlmetMin(const LeptonVector& l, const Met* met, float minVal)
{
  if( l.size() < 2 ) return false;
  const TLorentzVector *l0 = static_cast<TLorentzVector*>(l[0]);
  const TLorentzVector *l1 = static_cast<TLorentzVector*>(l[1]);
  float mtL0Met = SusyPlotter::transverseMass(*l0, met->lv());
  float mtL1Met = SusyPlotter::transverseMass(*l1, met->lv());
  return (mtL0Met < mtL1Met ? mtL0Met : mtL1Met) > minVal;
}
//----------------------------------------------------------
bool SusySelection::passMT2(const LeptonVector& leptons, const Met* met, float cut)
{
  if( leptons.size() < 2 ) return false;
  TLorentzVector metlv = met->lv();
  TLorentzVector l0    = *leptons.at(0);
  TLorentzVector l1    = *leptons.at(1);
  return (SusySelection::computeMt2(l0, l1, metlv) > cut);
}
//----------------------------------------------------------
bool SusySelection::passHtMin(const LeptonVector& leptons,
                              const JetVector &jets,
                              const Met* met,
                              float minVal)
{ // DG : SusyNtTools::Meff uses all leptons, all jets, and met; is this what we want?
  return (minVal < SusyNtTools::Meff(leptons, jets, met));
}
//----------------------------------------------------------
bool SusySelection::passNlepMin(const LeptonVector &leptons, size_t minVal)
{
  // return (leptons.size() >= minVal);
  // DG we should define a m_signalLeptons2L instead
  // similar to Anyes' implementation
  size_t nLep=0;
  for(size_t i=0;i<leptons.size(); ++i){
    if(const Susy::Lepton* l = leptons[i])
      if(l->isMu() && fabs(l->Eta())>2.4) // 2L muon trig, see Anyes' email 2013-08-02
        return false;
    nLep++;
  }
  return nLep>=minVal;
}
//----------------------------------------------------------
bool SusySelection::passNj(const JetVector& jets, int minNj, int maxNj)
{
  int nj(numberOfCLJets(jets));
  return (minNj < nj && nj <= maxNj
	  && numberOfCBJets(jets) < 1);
}
//----------------------------------------------------------
bool SusySelection::passZtautauVeto(cvl_t& l, cvj_t& j, const Met* m, float widthZpeak)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  if(!m)         return false;
  float mZ0(91.);
  return abs(SusyPlotter::mZTauTau(*l[0], *l[1], m->lv()) - mZ0) > widthZpeak;
}
//----------------------------------------------------------
bool SusySelection::passZllVeto(cvl_t& l, float mllLo, float mllHi)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  float mll((*l[0] + *l[1]).M());
  return (mll<mllLo || mllHi<mll);
}
//----------------------------------------------------------
void swap(float &a, float &b) { float c(a); a=b; b=c; };
bool SusySelection::pass2LepPt(cvl_t& l, float minPt0, float minPt1)
{
  if(l.size()<2) return false;
  float pt0(l[0]->Pt()), pt1(l[1]->Pt());
  if(pt0 < pt1) swap(pt0, pt1);
  return (pt0>minPt0 && pt1>minPt1);
}
//----------------------------------------------------------
bool SusySelection::passPtllMin(cvl_t& l, float minPt)
{
  if(l.size()<2) return false;
  return TLorentzVector(*l[0] + *l[1]).Pt() > minPt;
}
//----------------------------------------------------------
bool SusySelection::passPtTot(cvl_t& l, cvj_t& j, const Met* m, float maxPtTot)
{
  if(!m) return false;
  TLorentzVector mlv(m->lv()), ll, jj;
  if (l.size()>0 && l[0]) ll += *l[0];
  if (l.size()>1 && l[1]) ll += *l[1];
  if (j.size()>0 && j[0]) jj += *j[0];
  if (j.size()>1 && j[1]) jj += *j[1];
  return (ll+jj+mlv).Pt() < maxPtTot;
}
//----------------------------------------------------------
bool SusySelection::passMllMax(const LeptonVector& l, float maxMll)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  return TLorentzVector(*l[0] + *l[1]).M() < maxMll;
}
//----------------------------------------------------------
bool SusySelection::passMllMin(const LeptonVector& l, float minVal)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  return minVal < TLorentzVector(*l[0] + *l[1]).M();
}
//----------------------------------------------------------
bool SusySelection::passDrllMax(const LeptonVector& l, float maxDr)
{
  if(l.size()<2) return false;
  return  l[0]->DeltaR(*(static_cast<TLorentzVector*>(l[1]))) < maxDr;
}
//-----------------------------------------
// Check MC Lepton
//-----------------------------------------
bool SusySelection::isRealLepton(const Lepton* lep)
{
  // Updated way of handling real and fake leptons using LeptonTruthTools
  return (lep->truthType == RecoTruthMatch::PROMPT);
}
//-----------------------------------------
bool SusySelection::isFakeLepton(const Lepton* lep)
{
  return !isRealLepton(lep);
}
//-----------------------------------------
bool SusySelection::isConvLepton(const Lepton* lep)
{
  bool isConverted(lep->truthType == RecoTruthMatch::CONV);
  bool isChargeFlip(lep->isEle() ?
                    static_cast<const Electron*>(lep)->isChargeFlip : false);
  return isConverted && !isChargeFlip;
}
//-----------------------------------------
bool SusySelection::isHFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::HF);
}
//-----------------------------------------
bool SusySelection::isLFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::LF);
}
//-----------------------------------------
bool SusySelection::isTrueDilepton(const LeptonVector &leptons)
{
  // Maybe not 100% kosher, but I just want to make sure I am not
  // double counting, so I require dilepton events to be real
  if( leptons.size() != 2 ) return false;
  bool l0_real = isRealLepton(leptons[0]);
  bool l1_real = isRealLepton(leptons[1]);
  bool l0_cf = leptons[0]->isEle() ? ((Electron*) leptons[0])->isChargeFlip : false;
  bool l1_cf = leptons[1]->isEle() ? ((Electron*) leptons[1])->isChargeFlip : false;
  return l0_real && l1_real && !l0_cf && !l1_cf; // ignoring charge flip
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
bool SusySelection::passEleD0S(const LeptonVector &leptons, float maxVal)
{
  for(size_t i=0; i<leptons.size(); ++i){
      const Susy::Lepton* l = leptons[i];
      if(l->isEle() && (fabs(l->d0Sig(true)) > maxVal)) return false;
  } // end for(i)
  return true;
}

//-----------------------------------------
float SusySelection::getEvtWeight(const LeptonVector& leptons,
                                  bool includeBTag, bool includeTrig)
{
  float weight = 1.0;
  bool useSumwMap(true);
  if( !nt.evt()->isMC ) return weight;
  // lumi, xs, sumw, pileup
  weight = (m_useXsReader ?
            computeEventWeightXsFromReader(LUMI_A_L) :
            SusyNtAna::getEventWeight(LUMI_A_L, useSumwMap));
  weight *= getPythiaBbCcScaleFactor(nt.evt()->mcChannel, leptons);
  weight *= getLeptonEff2Lep(leptons);
  float trigW = (includeTrig ? getTriggerWeight2Lep(leptons) : 1.0);
  float bTag  = (includeBTag ? getBTagWeight(nt.evt()) : 1.0);
  //cout<<"bTag : "<<bTag<<endl;
  return weight * trigW * bTag;
}
//-----------------------------------------
float SusySelection::getPythiaBbCcScaleFactor(uint datasetId, const LeptonVector &leptons) const
{
  bool isBb2el(datasetId==129135), isBb2mu(datasetId==129136);
  bool isCc2el(datasetId==147667), isCc2mu(datasetId==147668);
  bool isPythiaBbCcSample(isBb2el || isBb2mu || isCc2el || isCc2mu);
  bool rigthLepton(isBb2mu || isCc2mu ?
                   leptons[0]->isMu()  && leptons[1]->isMu() :
                   leptons[0]->isEle() && leptons[1]->isEle());
  return (isPythiaBbCcSample ? (rigthLepton ? 0.706074 : 0.0) : 1.0);

}
//-----------------------------------------
float SusySelection::getBTagWeight(const Event* evt)
{
  JetVector tempJets;
  const JetVector& jets = m_signalJets2Lep; //m_baseJets
  for(uint ij=0; ij<jets.size(); ++ij){
    Jet* jet = jets.at(ij);
    if( !(jet->Pt() > 20 && fabs(jet->Eta()) < 2.5) ) continue;
    //if( fabs(jet->Eta()) > JET_ETA_CUT  ) continue;
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
float SusySelection::getLeptonEff2Lep(const LeptonVector &leptons) const
{
  assert(leptons.size()>1);
  return leptons[0]->effSF * leptons[1]->effSF;
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
std::string lineCountersPerEventType(const float cnt[ET_N][WT_N],
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
  string (&lcpet)(const float cnt[ET_N][WT_N], int weightType, int colWidth) = lineCountersPerEventType;
  for(int w=0; w<WT_N; ++w){
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
float SusySelection::computeMt2(const TLorentzVector &l0, const TLorentzVector &l1,
				const TLorentzVector &met)
{
  double pTMiss[3] = {0.0, met.Px(), met.Py()};
  double pA[3]     = {0.0, l0.Px(), l0.Py()};
  double pB[3]     = {0.0, l1.Px(), l1.Py()};
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pA,pB,pTMiss);
  mt2_event.set_mn(0); // LSP mass = 0 is Generic
  return mt2_event.get_mt2();
}
//-----------------------------------------
float SusySelection::computeChargeFlipProb(LeptonVector &leptons, Met &met,
                                           uint systematic, // DG todo
                                           bool update4mom)
{
  cvl_t &ls = leptons;
  if(ls.size()<2 || !ls[0] || !ls[1] || !m_chargeFlip) return 0.0;
  Lepton *l0(ls[0]), *l1(ls[1]);
  int pdg0(pdgIdFromLep(l0)), pdg1(pdgIdFromLep(l1));
  TLorentzVector smearedLv0(*l0), smearedLv1(*l1);
  TVector2 smearedMet(met.lv().Px(), met.lv().Py());
  int sys(NtSys_NOM==systematic ? 0 : 0);
  //(DGSys_BKGMETHOD_UP==systematic ? +1 : // DG todo : implement syst
  // (DGSys_BKGMETHOD_DN==systematic ? -1 : 0)));
  float flipProb(m_chargeFlip->OS2SS(pdg0, &smearedLv0, pdg1, &smearedLv1,
                                     &smearedMet, sys));
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
int SusySelection::pdgIdFromLep(const Lepton *l)
{
  // particles have positive codes, see doi:10.1146/annurev.ns.25.120175.003011
  int kPel(+11), kAel(-11), kPmu(+13), kAmu(-13), kUnknown(0);
  if     (l->isEle()) return (l->q < 0 ? kPel : kAel);
  else if(l->isMu() ) return (l->q < 0 ? kPmu : kAmu);
  else                return kUnknown;
}
//-----------------------------------------
void SusySelection::resetAllCounters()
{
  for(int w=0; w<WT_N; ++w){// Loop over weight types
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
  chargeFlipInput += "/../ChargeFlip/data/chargeFlip.root";
  m_chargeFlip = new chargeFlip(chargeFlipInput);
  if(m_dbg) m_chargeFlip->printSettings();
}
//-----------------------------------------
