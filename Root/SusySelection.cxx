#include <iomanip>
#include <cassert>
#include <cmath> // isnan
#include <iomanip> // setw
#include <sstream>
#include "TCanvas.h"
#include "SusyTest0/SusySelection.h"
#include "SusyTest0/SusyPlotter.h"

#include "Mt2/mt2_bisect.h"
#include "LeptonTruthTools/RecoTruthMatch.h" // provides RecoTruthMatch::
#include "SusyNtuple/WhTruthExtractor.h"

using namespace std;
using namespace Susy;

SusySelection::SusySelection() :
  m_susyObj(NULL),
  m_xsReader(NULL),
  m_trigObj(NULL),
  m_useMCTrig(false),
  m_fileName("default"),
  m_w(1.0),
  m_do1fb(false),
  m_doAD(false),
  m_useXsReader(false),
  m_xsFromReader(-1.0),
  m_dumpCounts(true),
  m_nLepMin(2),
  m_nLepMax(2),
  m_cutNBaseLep(true),
  m_ET(ET_Unknown)
{
  resetAllCounters();
  setAnaType(Ana_2Lep);
}
void SusySelection::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  if(m_dbg) cout << "SusySelection::Begin" << endl;

  string per = "Moriond";
  if(m_do1fb) per = "A-B3";
  if(m_doAD)  per = "A-D7";
  m_trigObj = new DilTrigLogic(per,false/*No Reweight Utils!*/);
  if(m_useMCTrig) m_trigObj->useMCTrigger();
  if( m_useXsReader ){
    m_xsReader = new XSReader();
    m_xsReader->setDebug(m_dbg);
    m_xsReader->LoadXSInfo();
  } // end if(m_useXsReader)
}

/*--------------------------------------------------------------------------------*/
// Main process loop function
/*--------------------------------------------------------------------------------*/
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
  selectObjects();
  if(!selectAnaEvent(m_signalLeptons, m_baseLeptons)) return kTRUE;
  bool count(true), includeLepSF(false);
  // Count SS and OS
  if(sameSign(m_signalLeptons))     increment(n_pass_ss[m_ET], includeLepSF);
  if(oppositeSign(m_signalLeptons)) increment(n_pass_os[m_ET], includeLepSF);
  // Check Signal regions
  passSR6(m_baseLeptons, m_signalJets2Lep, m_met, count);
  passSR7(m_baseLeptons, m_signalJets2Lep, m_met, count);
  passSR8(m_baseLeptons, m_signalJets2Lep, m_met, count);
  passSR9(m_baseLeptons, m_signalJets2Lep, m_met, count);
  return kTRUE;
}
void SusySelection::Terminate()
{
  SusyNtAna::Terminate();
  if(m_dbg) cout << "SusySelection::Terminate" << endl;
  if(m_dumpCounts)
    dumpEventCounters();
  if(m_xsReader) delete m_xsReader;
}

/*--------------------------------------------------------------------------------*/
// Full event selection
/*--------------------------------------------------------------------------------*/
bool SusySelection::selectEvent(bool doMll)
{
  if(m_dbg) cout << "SusySelection::selectEvent" << endl;
  // Basic event cuts

  int flag = nt.evt()->cutFlags[NtSys_NOM];
  int hdec = nt.evt()->hDecay;
  JetVector &jets = m_baseJets;
  const Susy::Met *met = m_met;
  uint run = nt.evt()->run;
  bool mc = nt.evt()->isMC;
  if(passGRL        (flag           ))  { increment(n_pass_Grl     ); } else { return false; }
  if(passLarErr     (flag           ))  { increment(n_pass_LarErr  ); } else { return false; }
  if(passTileErr    (flag           ))  { increment(n_pass_TileErr ); } else { return false; }
  if(passTTCVeto    (flag           ))  { increment(n_pass_TTCVeto ); } else { return false; }
  if(passGoodVtx    (flag           ))  { increment(n_pass_GoodVtx ); } else { return false; }
  if(passTileTripCut(flag           ))  { increment(n_pass_TileTrip); } else { return false; }
  if(passHfor       (               ))  {                           ; } else { return false; }
  if(passLAr        (flag           ))  { increment(n_pass_LAr     ); } else { return false; }
  if(passBadJet     (flag           ))  { increment(n_pass_BadJet  ); } else { return false; }
  if(passDeadRegions(jets,met,run,mc))  { increment(n_pass_FEBCut  ); } else { return false; }
  if(passBadMuon    (flag           ))  { increment(n_pass_BadMuon ); } else { return false; }
  if(passCosmic     (flag           ))  { increment(n_pass_Cosmic  ); } else { return false; }
  if(passHotSpot    (flag           ))  {                           ; } else { return false; }
  //if(passHtautauVeto(hdec)) { increment(n_pass_HttVeto ); } else { return false; }
  if( !nt.evt()->passMllForAlpgen ) return false;
  if(doMll && m_baseLeptons.size() == 2){
    if( Mll(m_baseLeptons[0], m_baseLeptons[1]) < 20 )
      return false;
  }
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::selectAnaEvent(const LeptonVector& leptons, const LeptonVector& baseLeps)
{

  if(m_dbg) cout << "SusySelection::selectAnaEvent" << endl;
  if( !selectEvent() )                     return false;
  // Lepton Analysis cuts
  if( baseLeps.size() < 2 )                return false;
  increment(n_pass_atleast2Lep);
  if(baseLeps.size() != 2)                 return false;
  //if(!(isFakeLepton(baseLeps[0]) || isFakeLepton(baseLeps[1]))) return false;
  increment(n_pass_exactly2Lep);

  m_ET = getDiLepEvtType(baseLeps);
  if(m_ET == ET_me) m_ET = ET_em;

  if( !passNLepCut(leptons) )              return false;

  // only check trigger for MC
  //if( !nt.evt()->isMC && !passTrigger(baseLeps) ) return false;
  if( !passTrigger(baseLeps) ) return false;

  if(m_ET == ET_mm){
    float trigW = m_trigObj->getTriggerWeight(m_baseLeptons,
					      nt.evt()->isMC,
					      m_met->Et,
					      m_signalJets2Lep.size(),
					      nt.evt()->nVtx,
					      NtSys_NOM);
    if(trigW < 0){
      cout<<"Run: "<<nt.evt()->run<<" Event: "<<nt.evt()->event<<" Weight: "<<trigW<<endl;
      m_baseLeptons[0]->print();
      m_baseLeptons[1]->print();
    }
  }
  return true;
}

/*--------------------------------------------------------------------------------*/
bool SusySelection::passSR6base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return oppositeSign(leptons) && sameFlavor(leptons);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passSR7base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return oppositeSign(leptons) && oppositeFlavor(leptons);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passSR8base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return sameSign(leptons) && sameFlavor(leptons);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passSR9base(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count)
{
  return sameSign(leptons) && oppositeFlavor(leptons);
}

/*--------------------------------------------------------------------------------*/
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
  if( passMETRel(met, leptons, jets) )     {if(count) increment(n_pass_SR6METRel      [m_ET], lepSf, bSf);}
  else return false;
  if( passMtLlmetMin(leptons, met) )       {if(count) increment(n_pass_SR6MtLlmetMin  [m_ET], lepSf, bSf);}
  else return false;
  if( passMtMinlmetMin(leptons, met) )     {if(count) increment(n_pass_SR6MtMinlmetMin[m_ET], lepSf, bSf);}
  else return false;
  if( passZtautauVeto(leptons, jets, met)) {if(count) increment(n_pass_SR6ZtautauVeto [m_ET], lepSf, bSf);}
  else return false;
  if( count )                               increment(n_pass_SR6[m_ET],lepSf, bSf);
  return true;
}
/*--------------------------------------------------------------------------------*/
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
  if( passMETRel(met, leptons, jets) )     {if(count) increment(n_pass_SR7METRel      [m_ET], lepSf, bSf);}
  else return false;
  if( passMtLlmetMin(leptons, met) )       {if(count) increment(n_pass_SR7MtLlmetMin  [m_ET], lepSf, bSf);}
  else return false;
  if( passMtMinlmetMin(leptons, met) )     {if(count) increment(n_pass_SR7MtMinlmetMin[m_ET], lepSf, bSf);}
  else return false;
  if( passZtautauVeto(leptons, jets, met)) {if(count) increment(n_pass_SR7ZtautauVeto [m_ET], lepSf, bSf);}
  else return false;
  //if( !passPtTot      (leptons, jets, met)) return false;
  if( count )                               increment(n_pass_SR7[m_ET],lepSf, bSf);
  return true;
}
/*--------------------------------------------------------------------------------*/
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
  if( passMETRel(met,leptons,jets,50.) ) {if( count ) increment(n_pass_SR8metr[m_ET],lepSf, bSf);}
  else return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
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
  if( passMETRel(met,leptons,jets,50.) ) {if( count ) increment(n_pass_SR9metr[m_ET],lepSf, bSf);}
  else return false;
  return true;
}

/*--------------------------------------------------------------------------------*/
// Generic cuts
/*--------------------------------------------------------------------------------*/
bool SusySelection::passHfor()
{
  if(nt.evt()->hfor == 4 ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passNLepCut(const LeptonVector& leptons)
{
  uint nLep = leptons.size();
  if(m_nLepMin>=0 && nLep < m_nLepMin) return false;
  if(m_nLepMax>=0 && nLep > m_nLepMax) return false;
  increment(n_pass_signalLep); //+=m_w;
  increment(n_pass_flavor[m_ET]);

  // To stay inline with Anders' cutflow
  if(m_ET == ET_me) increment(n_pass_flavor[ET_em],true);
  if(m_ET == ET_em) increment(n_pass_flavor[ET_me],true);

  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passNBaseLepCut(const LeptonVector& baseLeptons)
{
  if(m_cutNBaseLep){
    uint nLep = baseLeptons.size();
    if(m_nLepMin>=0 && nLep < m_nLepMin) return false;
    if(m_nLepMax>=0 && nLep > m_nLepMax) return false;
  }
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passTrigger(const LeptonVector& leptons)
{
  if(leptons.size() != 2) return false;
  bool passEvtTrig   = m_trigObj->passDilEvtTrig(leptons, m_met->Et, nt.evt());
  bool passTrigMatch = m_trigObj->passDilTrigMatch(leptons, m_met->Et, nt.evt());
  if( passEvtTrig ){ increment(n_pass_evtTrig[m_ET],true); }
  if( passEvtTrig && passTrigMatch){
    increment(n_pass_trigMatch[m_ET],true);
    return true;
  }
  return false;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::sameFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return (leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::oppositeFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return !(leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::sameSign(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return leptons.at(0)->q * leptons.at(1)->q > 0;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::oppositeSign(const LeptonVector& leptons)
{
  return !(sameSign(leptons));
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passMll(const LeptonVector& leptons, float mll)
{
  if( leptons.size() < 2 ) return false;
  if( (*leptons.at(0) + *leptons.at(1)).M() < mll ) return false;
  increment(n_pass_mll[m_ET]);
  return true;
}
bool SusySelection::passHtautauVeto(int hdecay)
{
  return (hdecay!=WhTruthExtractor::kPtauAtau);
}
/*--------------------------------------------------------------------------------*/
// Signal region cuts
/*--------------------------------------------------------------------------------*/
bool SusySelection::passJetVeto(const JetVector& jets)
{
  // Require no light, b, or forward jets
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 + N_B20 + N_F30 == 0);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passbJetVeto(const JetVector& jets)
{
  // Reject if there is a b jet using 2L definition
  int N_B20 = numberOfCBJets(jets);
  return (N_B20 == 0);
}
/*--------------------------------------------------------------------------------*/
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
/*--------------------------------------------------------------------------------*/
bool SusySelection::passeq2Jet(const JetVector& jets)
{
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L20 == 2 && N_B20 + N_F30 == 0);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passge2JetWoutFwVeto(const JetVector& jets)
{
  return (numberOfCLJets(jets) >= 2 && numberOfCBJets(jets) < 1);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passeq2JetWoutFwVeto(const JetVector& jets)
{
  return (numberOfCLJets(jets) == 2 && numberOfCBJets(jets) < 1);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passZVeto(const LeptonVector& leptons, float Zlow, float Zhigh)
{
  if( leptons.size() < 2 )   return false;
  //if( !sameFlavor(leptons) ) return true;
  float mll = (*leptons.at(0) + *leptons.at(1)).M();
  if( Zlow < mll && mll < Zhigh ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passMETRel(const Met *met, const LeptonVector& leptons,
				 const JetVector& jets, float metMax){
  return (getMetRel(met,leptons,jets) > metMax);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passdPhi(TLorentzVector v0, TLorentzVector v1, float cut)
{
  return v0.DeltaPhi(v1) > cut;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passMtLlmetMin(const LeptonVector& l, const Met* met, float minVal)
{
  if( l.size() < 2 ) return false;
  return SusyPlotter::transverseMass(TLorentzVector(*l[0]+*l[1]), met->lv()) > minVal;
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
  if(l.size()<2) return false;
  return TLorentzVector(*l[0] + *l[1]).M() < maxMll;
}
//----------------------------------------------------------
bool SusySelection::passDrllMax(const LeptonVector& l, float maxDr)
{
  if(l.size()<2) return false;
  return  l[0]->DeltaR(*(static_cast<TLorentzVector*>(l[1]))) < maxDr;
}
/*--------------------------------------------------------------------------------*/
// Check MC Lepton
/*--------------------------------------------------------------------------------*/
bool SusySelection::isRealLepton(const Lepton* lep)
{
  // Updated way of handling real and fake leptons using LeptonTruthTools
  return (lep->truthType == RecoTruthMatch::PROMPT);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isFakeLepton(const Lepton* lep)
{
  return !isRealLepton(lep);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isConvLepton(const Lepton* lep)
{
  bool isConv       = lep->truthType == RecoTruthMatch::CONV;
  bool isChargeFlip =  lep->isEle() ? static_cast<const Electron*>(lep)->isChargeFlip : false;
  return isConv && !isChargeFlip;
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isHFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::HF);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isLFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::LF);
}
/*--------------------------------------------------------------------------------*/
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
/*--------------------------------------------------------------------------------*/
float SusySelection::getEvtWeight(const LeptonVector& leptons, bool includeBTag, bool includeTrig)
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
  cout<<"bTag : "<<bTag<<endl;
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
/*--------------------------------------------------------------------------------*/
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
  return bTagSF(evt, tempJets, true, "MV1", "0_3511", MV1_80, BTag_NOM);
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
  return oss.str();
}
// helper function: for a given weight type, write line with counts for each event type
std::string lineCountersPerEventType(const float cnt[ET_N][WT_N], int weightType, int colWidth){
  std::ostringstream oss;
  for(int i=0; i<ET_N-1; ++i)
    oss<<std::setw(colWidth)<<cnt[i][weightType];
  return oss.str();
}
void SusySelection::dumpEventCounters()
{
  string v_ET[] = {"ee","mm","em","me"};
  string v_WT[] = {"Raw","Event","Pileup","Pileup A-B3",
                   "LeptonSF","btagSF","TrigSF","All A-B3", "All A-E"};
  int colWidth(10);
  int &cw = colWidth;
  int nCols(ET_N-1);
  string topRule(nCols*colWidth, '*');
  string midRule(nCols*colWidth, '-');
  // define a function reference to shorten lines
  string (&lcpet)(const float cnt[ET_N][WT_N], int weightType, int colWidth) = lineCountersPerEventType;
  for(int w=0; w<WT_N; ++w){
    cout<<topRule                                                         <<endl;
    cout<<"Event counts for weight: "<< v_WT             [w]              <<endl;
    cout<<midRule                                                         <<endl;
    cout<<"read in:              : "<<n_readin           [w]              <<endl;
    cout<<"pass GRL              : "<<n_pass_Grl         [w]              <<endl;
    cout<<"pass LarErr           : "<<n_pass_LarErr      [w]              <<endl;
    cout<<"pass TileErr          : "<<n_pass_TileErr     [w]              <<endl;
    cout<<"pass TTCVeto          : "<<n_pass_TTCVeto     [w]              <<endl;
    cout<<"pass GoodVtx          : "<<n_pass_GoodVtx     [w]              <<endl;
    cout<<"pass TileTripCut      : "<<n_pass_TileTrip    [w]              <<endl;
    cout<<"pass LAr:             : "<<n_pass_LAr         [w]              <<endl;
    cout<<"pass BadJet:          : "<<n_pass_BadJet      [w]              <<endl;
    cout<<"pass FEB:             : "<<n_pass_FEBCut      [w]              <<endl;
    cout<<"pass BadMu:           : "<<n_pass_BadMuon     [w]              <<endl;
    cout<<"pass Cosmic:          : "<<n_pass_Cosmic      [w]              <<endl;
    cout<<"pass Htautau veto     : "<<n_pass_HttVeto     [w]              <<endl;
    cout<<"   ------  Start Comparison Here ------ "                      <<endl;
    cout<<"pass atleast 2        : "<<n_pass_atleast2Lep [w]              <<endl;
    cout<<"pass exactly 2        : "<<n_pass_exactly2Lep [w]              <<endl;
    cout<<"pass nSigLep          : "<<n_pass_signalLep   [w]              <<endl;
    cout<<midRule                                                         <<endl;
    cout<<"For dilepton type     : "<<lineLabelsPerEventType(v_ET, cw)    <<endl;
    cout<<midRule                                                         <<endl;
    cout<<"pass flavor:          : "<<lcpet(n_pass_flavor         , w, cw)<<endl;
    cout<<"pass evt trig:        : "<<lcpet(n_pass_evtTrig        , w, cw)<<endl;
    cout<<"pass trig match:      : "<<lcpet(n_pass_trigMatch      , w, cw)<<endl;
    cout<<"pass OS:              : "<<lcpet(n_pass_os             , w, cw)<<endl;
    cout<<"pass SS:              : "<<lcpet(n_pass_ss             , w, cw)<<endl;
    cout<<midRule                                                         <<endl;
    cout<<"pass SR6 sign:        : "<<lcpet(n_pass_SR6sign        , w, cw)<<endl;
    cout<<"pass SR6 flavor:      : "<<lcpet(n_pass_SR6flav        , w, cw)<<endl;
    cout<<"pass SR6 >=1j:        : "<<lcpet(n_pass_SR6ge1j        , w, cw)<<endl;
    cout<<"pass SR6 >=2j:        : "<<lcpet(n_pass_SR6ge2j        , w, cw)<<endl;
    cout<<"pass SR6 DrllMax      : "<<lcpet(n_pass_SR6DrllMax     , w, cw)<<endl;
    cout<<"pass SR6 PtllMin      : "<<lcpet(n_pass_SR6PtllMin     , w, cw)<<endl;
    cout<<"pass SR6 MllMax       : "<<lcpet(n_pass_SR6MllMax      , w, cw)<<endl;
    cout<<"pass SR6 METRel       : "<<lcpet(n_pass_SR6METRel      , w, cw)<<endl;
    cout<<"pass SR6 MtLlmetMin   : "<<lcpet(n_pass_SR6MtLlmetMin  , w, cw)<<endl;
    cout<<"pass SR6 MtMinlmetMin : "<<lcpet(n_pass_SR6MtMinlmetMin, w, cw)<<endl;
    cout<<"pass SR6 ZtautauVeto  : "<<lcpet(n_pass_SR6ZtautauVeto , w, cw)<<endl;
    cout<<"pass SR6              : "<<lcpet(n_pass_SR6            , w, cw)<<endl;
    cout<<midRule                                                         <<endl;
    cout<<"pass SR7 sign:        : "<<lcpet(n_pass_SR7sign        , w, cw)<<endl;
    cout<<"pass SR7 flavor:      : "<<lcpet(n_pass_SR7flav        , w, cw)<<endl;
    cout<<"pass SR7 >=1j:        : "<<lcpet(n_pass_SR7ge1j        , w, cw)<<endl;
    cout<<"pass SR7 >=2j:        : "<<lcpet(n_pass_SR7ge2j        , w, cw)<<endl;
    cout<<"pass SR7 DrllMax      : "<<lcpet(n_pass_SR7DrllMax     , w, cw)<<endl;
    cout<<"pass SR7 PtllMin      : "<<lcpet(n_pass_SR7PtllMin     , w, cw)<<endl;
    cout<<"pass SR7 MllMax       : "<<lcpet(n_pass_SR7MllMax      , w, cw)<<endl;
    cout<<"pass SR7 METRel       : "<<lcpet(n_pass_SR7METRel      , w, cw)<<endl;
    cout<<"pass SR7 MtLlmetMin   : "<<lcpet(n_pass_SR7MtLlmetMin  , w, cw)<<endl;
    cout<<"pass SR7 MtMinlmetMin : "<<lcpet(n_pass_SR7MtMinlmetMin, w, cw)<<endl;
    cout<<"pass SR7 ZtautauVeto  : "<<lcpet(n_pass_SR7ZtautauVeto , w, cw)<<endl;
    cout<<"pass SR7              : "<<lcpet(n_pass_SR7            , w, cw)<<endl;
    cout<<midRule                                                         <<endl;
    cout<<"pass SR8 sign         : "<<lcpet(n_pass_SR8sign        , w, cw)<<endl;
    cout<<"pass SR8 flavor       : "<<lcpet(n_pass_SR8flav        , w, cw)<<endl;
    cout<<"pass SR8 >=1j         : "<<lcpet(n_pass_SR8ge1j        , w, cw)<<endl;
    cout<<"pass SR8 >=2j         : "<<lcpet(n_pass_SR8ge2j        , w, cw)<<endl;
    cout<<"pass SR8 METRel > 50  : "<<lcpet(n_pass_SR8metr        , w, cw)<<endl;
    cout<<midRule                                                         <<endl;
    cout<<"pass SR9 sign         : "<<lcpet(n_pass_SR9sign        , w, cw)<<endl;
    cout<<"pass SR9 flavor       : "<<lcpet(n_pass_SR9flav        , w, cw)<<endl;
    cout<<"pass SR9 >=1j         : "<<lcpet(n_pass_SR9ge1j        , w, cw)<<endl;
    cout<<"pass SR9 >=2j         : "<<lcpet(n_pass_SR9ge2j        , w, cw)<<endl;
    cout<<"pass SR9 METRel > 50  : "<<lcpet(n_pass_SR9metr        , w, cw)<<endl;
    cout<<midRule                                                         <<endl;
  }// end for(w)
}

/*--------------------------------------------------------------------------------*/
void SusySelection::dumpPreObjects()
{
  ElectronVector preElectrons = getPreElectrons(&nt, NtSys_NOM);
  MuonVector preMuons = getPreMuons(&nt, NtSys_NOM);
  JetVector preJets = getPreJets(&nt, NtSys_NOM);
  cout<<"Pre Electrons: "<<preElectrons.size()<<endl;
  for(uint ie=0; ie<preElectrons.size(); ++ie){ preElectrons[ie]->print(); }
  cout<<endl;
  cout<<"Pre Muons: "<<preMuons.size()<<endl;
  for(uint im=0; im<preMuons.size(); ++im){ preMuons[im]->print(); }
  cout<<endl;
  cout<<"Pre Jets: "<<preJets.size()<<endl;
  for(uint ij=0; ij<preJets.size(); ++ij){ preJets[ij]->print(); }
}
/*--------------------------------------------------------------------------------*/
// Dump Jets with more information
/*--------------------------------------------------------------------------------*/
void SusySelection::dumpJets()
{
  cout<<"Jets:"<<endl;
  for(uint j=0; j<m_signalJets2Lep.size(); ++j){ m_signalJets2Lep.at(j)->print(); }
}
/*--------------------------------------------------------------------------------*/
// Method for checking systematics
/*--------------------------------------------------------------------------------*/
void SusySelection::checkSys()
{

  bool electronSys = false;
  bool muonSys = false; //false;
  bool jetSys = false;
  bool metSys = true; //true;
  bool effSys = false; //true; //true;
  cout<<endl;
  cout<<"*********************************************"<<endl;
  cout<<"Run: "<<nt.evt()->run<<" Event: "<<nt.evt()->event<<endl;
  // Electron systematics
  if(electronSys){
    for(uint i=0; i<m_signalLeptons.size(); ++i){
      if( !m_signalLeptons.at(i)->isEle() ) continue;
      const Electron* lep = (Electron*) m_signalLeptons.at(i);
      cout<<"Electron Pt: "<<lep->Pt()<<" Eta: "<<lep->Eta()<<" Phi: "<<lep->Phi()<<endl;
      cout<<"\tEnergy Scale Z Up:            "<<lep->Pt() * lep->ees_z_up<<endl;
      cout<<"\tEnergy Scale Z Down:          "<<lep->Pt() * lep->ees_z_dn<<endl;
      cout<<"\tEnergy Scale Material Up:     "<<lep->Pt() * lep->ees_mat_up<<endl;
      cout<<"\tEnergy Scale Material Down:   "<<lep->Pt() * lep->ees_mat_dn<<endl;
      cout<<"\tEnergy Scale Presampler Up:   "<<lep->Pt() * lep->ees_ps_up<<endl;
      cout<<"\tEnergy Scale Presampler Down: "<<lep->Pt() * lep->ees_ps_dn<<endl;
      cout<<"\tEnergy Scale Low Pt Up:       "<<lep->Pt() * lep->ees_low_up<<endl;
      cout<<"\tEnergy Scale Low Pt Down:     "<<lep->Pt() * lep->ees_low_dn<<endl;
      cout<<"\tEnergy Resolution Up:         "<<lep->Pt() * lep->eer_up<<endl;
      cout<<"\tEnergy Resolution Down:       "<<lep->Pt() * lep->eer_dn<<endl;
      cout<<endl;
    }// end loop over leptons
  } // end electron sys
  // Muon Systematics
  if(muonSys){
    for(uint i=0; i<m_signalLeptons.size(); ++i){
      if( m_signalLeptons.at(i)->isEle() ) continue;
      const Muon* lep = (Muon*) m_signalLeptons.at(i);
      cout<<"Muon Pt: "<<lep->Pt()<<" Eta: "<<lep->Eta()<<" Phi: "<<lep->Phi()<<endl;
      cout<<"\tMuon MS Up:   "<<lep->Pt() * lep->ms_up<<endl;
      cout<<"\tMuon MS Down: "<<lep->Pt() * lep->ms_dn<<endl;
      cout<<"\tMuon ID Up:   "<<lep->Pt() * lep->id_up<<endl;
      cout<<"\tMuon ID Down: "<<lep->Pt() * lep->id_dn<<endl;
      cout<<endl;
    }
  }// end muon sys

  // Jet Systematics
  if(jetSys){
    if(m_signalJets2Lep.size() > 0) cout<<"-----------------------------------"<<endl;
    for(uint i=0; i<m_signalJets2Lep.size(); ++i){
      const Jet* jet = m_signalJets2Lep.at(i);
      cout<<"Jet Pt: "<<jet->Pt()<<" Eta: "<<jet->Eta()<<" Phi: "<<jet->Phi()<<endl;
      cout<<"\t JES Up:   "<<jet->Pt() * jet->jes_up<<endl;
      cout<<"\t JES Down: "<<jet->Pt() * jet->jes_dn<<endl;
      cout<<"\t JER:      "<<jet->Pt() * jet->jer<<endl;
      cout<<endl;
    }// end loop over jets
  }// end if jet sys

  // Print Met info for all the shifts
  if(metSys){
    #define printMet(met) cout<<" Et: "<<met->Et<<" Phi: "<<met->phi<<endl;
    cout<<"Met Systematics:"<<endl;
    for(int i=NtSys_SCALEST_UP; i<NtSys_TRIGSF_EL_UP; ++i){
      cout<<"\t"<<i<<endl;
      const Met* tmp_met = getMet(&nt, (SusyNtSys) i);
      if(i == NtSys_NOM)            { cout<<"\tNominal:                         "; printMet(tmp_met); }
      else if(i == NtSys_EES_Z_UP)  { cout<<"\tEl Energy Scale Z UP:            "; printMet(tmp_met); }
      else if(i == NtSys_EES_Z_DN)  { cout<<"\tEl Energy Scale Z Down:          "; printMet(tmp_met); }
      else if(i == NtSys_EES_MAT_UP){ cout<<"\tEl Energy Scale Material UP:     "; printMet(tmp_met); }
      else if(i == NtSys_EES_MAT_DN){ cout<<"\tEl Energy Scale Material Down:   "; printMet(tmp_met); }
      else if(i == NtSys_EES_PS_UP) { cout<<"\tEl Energy Scale PreSampler UP:   "; printMet(tmp_met); }
      else if(i == NtSys_EES_PS_DN) { cout<<"\tEl Energy Scale PreSampler Down: "; printMet(tmp_met); }
      else if(i == NtSys_EES_LOW_UP){ cout<<"\tEl Energy Scale Low Pt UP:       "; printMet(tmp_met); }
      else if(i == NtSys_EES_LOW_DN){ cout<<"\tEl Energy Scale Low Pt Down:     "; printMet(tmp_met); }
      else if(i == NtSys_EER_UP)    { cout<<"\tEl Energy Resolution UP:         "; printMet(tmp_met); }
      else if(i == NtSys_EER_DN)    { cout<<"\tEl Energy Resolution Down:       "; printMet(tmp_met); }
      else if(i == NtSys_JES_UP)    { cout<<"\tJet Energy Scale UP:             "; printMet(tmp_met); }
      else if(i == NtSys_JES_DN)    { cout<<"\tJet Energy Scale Down:           "; printMet(tmp_met); }
      else if(i == NtSys_JER)       { cout<<"\tJet Energy Resolution:           "; printMet(tmp_met); }
      else if(i == NtSys_SCALEST_UP){ cout<<"\tMet Scale Soft Term UP:          "; printMet(tmp_met); }
      else if(i == NtSys_SCALEST_DN){ cout<<"\tMet Scale Soft Term Down:        "; printMet(tmp_met); }
      else if(i == NtSys_RESOST)    { cout<<"\tMet Resolution Soft Term:        "; printMet(tmp_met); }

    }
    #undef printMet
  }// end if met sys

  if(effSys){
    float elEff = 1;
    float elEffUp = 1;
    float elEffDn = 1;
    float muEff = 1;
    float muEffUp = 1;
    float muEffDn = 1;
    for(uint i=0; i<m_signalLeptons.size(); ++i){
      const Lepton* lep = m_signalLeptons.at(i);
      if(lep->isEle()){
        elEff   *= lep->effSF;
        elEffUp *= (lep->effSF + lep->errEffSF);
        elEffDn *= (lep->effSF - lep->errEffSF);
      }
      else{
        muEff   *= lep->effSF;
        muEffUp *= (lep->effSF + lep->errEffSF);
        muEffDn *= (lep->effSF - lep->errEffSF);
        cout<<"Error: "<<lep->errEffSF<<endl;
      }
    }// end loop over leptons

    cout<<"Efficiency Systematics:"<<endl;
    cout<<"\tNominal:     "<<elEff*muEff<<endl;
    cout<<"\tEl Eff Up:   "<<elEffUp * muEff<<endl;
    cout<<"\tEl Eff Down: "<<elEffDn * muEff<<endl;
    cout<<"\tMu Eff Up:   "<<elEff   * muEffUp<<endl;
    cout<<"\tMu Eff Down: "<<elEff   * muEffDn<<endl;
  }// end if efficiency sys
}
/*--------------------------------------------------------------------------------*/
// Debug event
/*--------------------------------------------------------------------------------*/
bool SusySelection::debugEvent()
{
  uint run = nt.evt()->run;
  uint evt = nt.evt()->event;
  if(run == 0 && evt == 0) return true;
  return false;
}
/*--------------------------------------------------------------------------------*/
// Dump interesting events
/*--------------------------------------------------------------------------------*/
void SusySelection::dumpInterestingEvents(const LeptonVector& leptons,
                                          const JetVector& jets,
                                          const Met* met)
{
  // Look for 1 jet, SS muons, 90-120 GeV
  if(!sameSign(leptons)) return;
  if(!(leptons[0]->isMu() && leptons[1]->isMu())) return;
  float mll = Mll(leptons[0],leptons[1]);
  if(!(100 < mll && mll < 110)) return;
  if(jets.size() != 1) return;
  if( met->Et < 40 ) return;
  out<<"-----------------------------------------------"<<endl;
  out<<"Run: "<<nt.evt()->run<<" Event: "<<nt.evt()->event<<endl;
  out<<"Mll "<<mll<<" Met: "<<met->Et<<endl;
  out<<"Signal Leptons: "<<endl;
  for(uint i=0; i<leptons.size(); ++i){
    out<<"Lepton "<<i<<endl;
    printLep(leptons[i]);
  }
  out<<"++++++++++++++++++++++"<<endl;
  out<<"Pre Leptons: "<<endl;
  ElectronVector elecs = getPreElectrons(&nt, NtSys_NOM);
  MuonVector muons     = getPreMuons(&nt, NtSys_NOM);
  TauVector taus       = getPreTaus(&nt, NtSys_NOM);
  for(uint i=0; i<elecs.size(); ++i){ out<<"Pre Electron: "<<i<<endl; printLep((Lepton*) elecs[i]); }
  for(uint i=0; i<muons.size(); ++i){ out<<"Pre Muons: "<<i<<endl; printLep((Lepton*) muons[i]); }
  for(uint i=0; i<taus.size(); ++i){ out<<"Pre Taus: "<<i<<endl; printLep((Lepton*) taus[i]); }
  out<<"++++++++++++++++++++++"<<endl;
  out<<"N Signal Jets"<<jets.size()<<endl;
  for(uint i=0; i<jets.size(); ++i){ out<<"Jet: "<<i<<endl; printJet(jets[i]); }
  out<<"++++++++++++++++++++++"<<endl;
  JetVector prejets = getPreJets(&nt, NtSys_NOM);
  out<<"N Pre Jets"<<prejets.size()<<endl;
  for(uint i=0; i<prejets.size(); ++i){ out<<"Jet: "<<i<<endl; printJet(prejets[i]); }
}
/*--------------------------------------------------------------------------------*/
void SusySelection::dumpTrigFlag(uint flag)
{
  out << "\tEF_e7_medium1               " << (flag & TRIG_e7_medium1)               << endl;
  out << "\tEF_e12Tvh_medium1           " << (flag & TRIG_e12Tvh_medium1)           << endl;
  out << "\tEF_e24vh_medium1            " << (flag & TRIG_e24vh_medium1)            << endl;
  out << "\tEF_e24vhi_medium1           " << (flag & TRIG_e24vhi_medium1)           << endl;
  out << "\tEF_2e12Tvh_loose1           " << (flag & TRIG_2e12Tvh_loose1)           << endl;
  out << "\tEF_e24vh_medium1_e7_medium1 " << (flag & TRIG_e24vh_medium1_e7_medium1) << endl;
  out << "\tEF_mu8                      " << (flag & TRIG_mu8)                      << endl;
  out << "\tEF_mu18_tight               " << (flag & TRIG_mu18_tight)               << endl;
  out << "\tEF_mu24i_tight              " << (flag & TRIG_mu24i_tight)              << endl;
  out << "\tEF_2mu13                    " << (flag & TRIG_2mu13)                    << endl;
  out << "\tEF_mu18_tight_mu8_EFFS      " << (flag & TRIG_mu18_tight_mu8_EFFS)      << endl;
  out << "\tEF_e12Tvh_medium1_mu8       " << (flag & TRIG_e12Tvh_medium1_mu8)       << endl;
  out << "\tEF_mu18_tight_e7_medium1    " << (flag & TRIG_mu18_tight_e7_medium1)    << endl;
}
/*--------------------------------------------------------------------------------*/
void SusySelection::printLep(const Lepton* lep)
 {
  out<<"\tLepton is Muon: "<<lep->isMu()<<" Charge: "<<lep->q<<endl;
  if(lep->isMu())out<<"\tCombined: "<<((Muon*)lep)->isCombined<<endl;
  out<<"\tPt: "<<lep->Pt()<<" Eta: "<<lep->Eta()<<" Phi: "<<lep->Phi()<<endl;
  out<<"\tptcone20: "<<lep->ptcone20<<" ptcone30: "<<lep->ptcone30<<endl;
  out<<"\td0: "<<lep->d0<<" d0err: "<<lep->errD0<<endl;
}
/*--------------------------------------------------------------------------------*/
void SusySelection::printJet(const Jet* jet)
{
  out<<"\tPt: "<<jet->Pt()<<" Eta: "<<jet->Eta()<<" Phi: "<<jet->Phi()<<endl;
  out<<"\tjvf: "<<jet->jvf<<" sv0: "<<jet->sv0
     <<" combNN: "<<jet->combNN<<" mv1: "<<jet->mv1<<endl;
}
/*--------------------------------------------------------------------------------*/
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
/*--------------------------------------------------------------------------------*/
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
    n_pass_HttVeto    [w] = 0;
    n_pass_atleast2Lep[w] = 0;
    n_pass_exactly2Lep[w] = 0;
    n_pass_signalLep  [w] = 0;
    for(int i=0; i<ET_N; ++i){ // loop over weight x channel.
      n_pass_flavor[i][w]   = 0;
      n_pass_evtTrig[i][w]     = 0;
      n_pass_trigMatch[i][w]   = 0;
      n_pass_mll[i][w]         = 0;
      n_pass_ss[i][w]          = 0;
      n_pass_os[i][w]          = 0;
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
    } // end for(i)
  } // end for(w)
}
