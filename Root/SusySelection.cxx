#include <iomanip>
#include <cassert>
#include "TCanvas.h"
#include "SusyTest0/SusySelection.h"
#include "SusyTest0/SusyPlotter.h"

#include "Mt2/mt2_bisect.h"
#include "LeptonTruthTools/RecoTruthMatch.h" // provides RecoTruthMatch::

using namespace std;
using namespace Susy;

/*--------------------------------------------------------------------------------*/
// SusySelection Constructor
/*--------------------------------------------------------------------------------*/
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
  // Loop over weight types
  for(int w=0; w<WT_N; ++w){
    n_readin[w]       = 0;
    n_pass_LAr[w]     = 0;
    n_pass_BadJet[w]  = 0;
    n_pass_BadMuon[w] = 0;
    n_pass_Cosmic[w]  = 0;
    n_pass_atleast2Lep[w] = 0;
    n_pass_exactly2Lep[w] = 0;
    n_pass_signalLep[w]   = 0;

    // The rest are channel specific.
    for(int i=0; i<ET_N; ++i){
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

      n_pass_SR7sign[i][w] = n_pass_SR7flav[i][w] = n_pass_SR7metr[i][w] = 0;
      n_pass_SR7ge1j[i][w] = n_pass_SR7ge2j[i][w] = n_pass_SR7eq2j[i][w] = 0;
      n_pass_SR7eq2jNfv[i][w] = n_pass_SR7ge2jNfv[i][w] = n_pass_SR7[i][w] = 0;

      n_pass_SR8sign[i][w] = n_pass_SR8flav[i][w] = n_pass_SR8metr[i][w] = 0;
      n_pass_SR8ge1j[i][w] = n_pass_SR8ge2j[i][w] = n_pass_SR8eq2j[i][w] = 0;
      n_pass_SR8eq2jNfv[i][w] = n_pass_SR8ge2jNfv[i][w] = n_pass_SR8[i][w] = 0;

      n_pass_SR9sign[i][w] = n_pass_SR9flav[i][w] = n_pass_SR9metr[i][w] = 0;
      n_pass_SR9ge1j[i][w] = n_pass_SR9ge2j[i][w] = n_pass_SR9eq2j[i][w] = 0;
      n_pass_SR9eq2jNfv[i][w] = n_pass_SR9ge2jNfv[i][w] = n_pass_SR9[i][w] = 0;
    }
  }// end loop over weight types

  //out.open("event.dump");

  //setAnaType(Ana_2Lep);
  //out.open("InterestingEventsWPreTaus.txt");
  //out.open("InterestingEventsNoJetReq.txt");
  out.open("InterestingEvents3Lep.txt");
  //out.open("SSMMInclusive.txt");
  //out.open("dump.txt");
}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
/*--------------------------------------------------------------------------------*/
void SusySelection::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  if(m_dbg) cout << "SusySelection::Begin" << endl;

  string per = "HCP";
  if(m_do1fb) per = "A-B3";
  if(m_doAD)  per = "A-D7";
  m_trigObj = new DilTrigLogic(per);
  if(m_useMCTrig) m_trigObj->useMCTrigger();

  if( m_useXsReader ){
    m_xsReader = new XSReader();
    m_xsReader->setDebug(m_dbg);
    m_xsReader->LoadXSInfo();
  } // end if(m_useXsReader)
  //setSelectTaus(true);
}

/*--------------------------------------------------------------------------------*/
// Main process loop function
/*--------------------------------------------------------------------------------*/
Bool_t SusySelection::Process(Long64_t entry)
{
  // Communicate tree entry number to SusyNtObject
  GetEntry(entry);
  clearObjects();
  m_ET = ET_Unknown;
  increment(n_readin);

  // Chain entry not the same as tree entry
  //static Long64_t chainEntry = -1;
  m_chainEntry++;
  if(m_dbg || m_chainEntry%50000==0)
  {
    cout << "**** Processing entry " << setw(6) << m_chainEntry
         << " run " << setw(6) << nt.evt()->run
         << " event " << setw(7) << nt.evt()->event << " ****" << endl;
  }

  // select signal objects
  selectObjects();
  //dumpBaselineObjects();
  // Check Event
  if(!selectAnaEvent(m_signalLeptons, m_baseLeptons)) return kTRUE;


  //--- NEW ---//
  // Dump intersting events
  //dumpInterestingEvents(m_signalLeptons,m_signalJets2Lep, m_met);
  //return kTRUE;


  // Count SS and OS
  if(sameSign(m_signalLeptons))     increment(n_pass_ss[m_ET], true);
  if(oppositeSign(m_signalLeptons)) increment(n_pass_os[m_ET], true);

  // Check Signal regions
  bool count(true);
  passSR7(m_baseLeptons, m_signalJets2Lep, m_met, count);

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called
/*--------------------------------------------------------------------------------*/
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
  if( !passHfor() )                 return false;
  if( !passLAr(flag) )              return false;
  increment(n_pass_LAr);
  if( !passBadJet(flag) )           return false;
  increment(n_pass_BadJet);
  if( !passBadMuon(flag) )          return false;
  increment(n_pass_BadMuon);
  if( !passCosmic(flag) )           return false;
  increment(n_pass_Cosmic);
  if( !passHotSpot(flag) )          return false;
  //--- NEW ---//
  //if( hasJetInBadFCAL(m_baseJets) ) return false;

  //--- NEW ---//
  if( !nt.evt()->passMllForAlpgen ) return false;
  if(doMll && m_baseLeptons.size() == 2){
    //cout<<"Checking mll..."<<endl;
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
    float trigW = m_trigObj->getTriggerWeight(m_baseLeptons,nt.evt()->isMC,NtSys_NOM);
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
  if( !passDrllMax(leptons) )               return false;
  if( !passPtllMin(leptons) )               return false;
  if( !passMllMax(leptons) )                return false;
  if( !passMETRelMin(met, leptons, jets) )  return false;
  if( !passMtLlmetMin(leptons, met) )       return false;
  if( !passMtMinlmetMin(leptons, met) )     return false;
  if( !passZtautauVeto(leptons, jets, met)) return false;
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
  // DG Mar 2013 : might want to add counters here as well
  if( !passDrllMax(leptons) )               return false;
  if( !passPtllMin(leptons) )               return false;
  if( !passMllMax(leptons) )                return false;
  if( !passMETRelMin(met, leptons, jets) )  return false;
  if( !passMtLlmetMin(leptons, met) )       return false;
  if( !passMtMinlmetMin(leptons, met) )     return false;
  if( !passZtautauVeto(leptons, jets, met)) return false;
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

  //cout<<"NLeptons: "<<leptons.size()<<endl;
  //for(uint i=0; i<leptons.size(); ++i)
  //cout<<"\tElectron: "<<leptons.at(i)->isEle()<<" "<<leptons.at(i)->trigFlags<<endl;

  if(leptons.size() != 2)
    return false;

  bool passEvtTrig   = m_trigObj->passDilEvtTrig(leptons, nt.evt());
  bool passTrigMatch = m_trigObj->passDilTrigMatch(leptons, nt.evt());

  if( passEvtTrig ){
    increment(n_pass_evtTrig[m_ET],true);
    // To stay inline with Anders' cutflow
    //if(m_ET == ET_me) increment(n_pass_flavor[ET_em]);
    //if(m_ET == ET_em) increment(n_pass_flavor[ET_me]);
  }
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
  //return (leptons.at(0)->isEle() && leptons.at(1)->isEle());
  //return (leptons.at(0)->isMu() && leptons.at(1)->isMu());
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

  /*
  bool failjet = false;
  for(uint i=0; i<jets.size(); ++i){
    const Jet* jet = jets.at(i);
    if( jet->Pt() < 30 ) continue;
    failjet = true;
    break;
  }

  if( failjet ) return false;
  return true;
  */
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::passbJetVeto(const JetVector& jets)
{

  // Reject if there is a b jet using 2L definition
  int N_B20 = numberOfCBJets(jets);
  return (N_B20 == 0);

  /*
  bool hasbjet = false;
  for(uint i=0; i<jets.size(); ++i){
    const Jet* jet = jets.at(i);
    if( jet->Pt() < 30 )       continue;
    if( !isBJet(jet, MV1_85) ) continue;
    hasbjet = true;
    break;
  }

  if( hasbjet ) return false;
  return true;
  */
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

  // N_L25 >=2 N_B20 = N_F30 = 0
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);

  return (N_L20 >=2 && N_B20 + N_F30 == 0);

  /*
  // Count jets above 30 GeV
  int njet = 0;
  for(uint j=0; j<jets.size(); ++j)
    if( jets.at(j)->Pt() > 30 )
      njet++;

  return (njet >= 2);
  */
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
  // N_L25 >=2 N_B20 = N_F30 = 0
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


  if( getMetRel(met,leptons,jets) < metMax ) return false;
  return true;
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

  // necessary variables
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
  return (lep->truthMatchType == RecoTruthMatch::PROMPT);

  // Code taken from Steve.  There seems to be an issue with Sherpa samples, so
  // need to handle those separately. Also just for clarification:
  // * mcOrigin = 9 -- Tau Lepton
  // * mcType   = 1 -- Unknown Electron
  // * mcType   = 2 -- Iso Electron
  // * mcType   = 5 -- Unknown Muon
  // * mcType   = 6 -- Iso Muon

  // Cut is sample dependent due to Sherpa classifications broken
  uint mcId = nt.evt()->mcChannel;

  // All tau leptons are classified as non-iso
  // I'm not sure why, yet, but for now I will treat them as real leptons.
  if(lep->mcOrigin == 9) return true;

  const int mcType = lep->mcType;

  // Sherpa diboson, assume all unknowns are real leptons
  // This is an approximation, but probably ok.
  if( (mcId>=126892 && mcId<=126895) || (mcId>=147770 && mcId<=147772) ||
      (mcId>=147774 && mcId<=147776)){
    if(lep->isEle()) return mcType == 1 || mcType == 2;
    else             return mcType == 5 || mcType == 6;
  }
  else{

    // 2-lep classifies everything as real if it
    // is from W, Z, tau, or top..
    //uint origin = lep->mcOrigin;
    //return origin == 9 || origin == 12 || origin == 13 || origin == 10;

    if(lep->isEle()) return mcType == 2;
    else             return mcType == 6;
  }


}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isFakeLepton(const Lepton* lep)
{
  return !isRealLepton(lep);
}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isConvLepton(const Lepton* lep)
{
  //return lep->mcOrigin == 5;
  bool isConv       = lep->truthMatchType == RecoTruthMatch::CONV;
  //bool isConv       = lep->mcOrigin == 5;
  bool isChargeFlip =  lep->isEle() ? ((Electron*) lep)->isChargeFlip : false;
  return isConv && !isChargeFlip;

}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isHFLepton(const Lepton* lep)
{

  return (lep->truthMatchType == RecoTruthMatch::HF);

  //uint origin = lep->mcOrigin;
  //return origin == 25 ||origin == 26 || origin == 27 || origin == 28 ||
  //origin == 29 || origin == 32 || origin == 33;



}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isLFLepton(const Lepton* lep)
{

  return (lep->truthMatchType == RecoTruthMatch::LF);

  // Steve's way:
  //bool isChargeFlip = lep->isEle() ? ((Electron*) lep)->isChargeFlip : false;
  //return isFakeLepton(lep) && !isConvLepton(lep) && !isHFLepton(lep) && !isChargeFlip;

  // 2-lep way:
  //uint origin = lep->mcOrigin;
  //return origin == 0 || origin == 23 || origin == 24 || origin == 30 ||
  //origin == 31 || origin == 34 || origin == 35;


}
/*--------------------------------------------------------------------------------*/
bool SusySelection::isTrueDilepton(const LeptonVector &leptons)
{

  //if( leptons.size() !=2 ) return false;
  //  bool l0real = leptons[0]->truthMatchType == RecoTruthMatch::PROMPT;
  //bool l1real = leptons[1]->truthMatchType == RecoTruthMatch::PROMPT;
  //return l0real && l1real;


  // Maybe not 100% kosher, but I just want to make sure I am not
  // double counting, so I require dilepton events to be real
  if( leptons.size() != 2 ) return false;

  bool l0_real = isRealLepton(leptons[0]);
  bool l1_real = isRealLepton(leptons[1]);
  bool l0_cf = leptons[0]->isEle() ? ((Electron*) leptons[0])->isChargeFlip : false;
  bool l1_cf = leptons[1]->isEle() ? ((Electron*) leptons[1])->isChargeFlip : false;

  return l0_real && l1_real && !l0_cf && !l1_cf; // ignoring charge flip
  //return (l0_real || l0_cf) && (l1_real || l1_cf);

}

/*--------------------------------------------------------------------------------*/
// Get Event weight
/*--------------------------------------------------------------------------------*/
float SusySelection::getEvtWeight(const LeptonVector& leptons, bool includeBTag, bool includeTrig,
				  bool doMediumpp)
{
  if( !nt.evt()->isMC ) return 1.;
  uint nl = leptons.size();
  float weight = 1;

  // lumi, xs, sumw, pileup
  if(m_do1fb) weight = getEventWeightAB3();
  else if(m_doAD)  weight = getEventWeight(LUMI_A_D);
  else weight = (m_useXsReader ? computeEventWeightXsFromReader(LUMI_A_E) : getEventWeight(LUMI_A_E));
  //if(m_do1fb) weight = getEventWeightFixed(nt.evt()->mcChannel, LUMI_A_B3);
  //else if(m_doAD)  weight = getEventWeightFixed(nt.evt()->mcChannel,LUMI_A_D);
  //else weight = getEventWeightFixed(nt.evt()->mcChannel,LUMI_A_E);

  // bbbar/ccbar scale factor
  uint chNum = nt.evt()->mcChannel;
  if(chNum == 129136 || chNum == 147668){
    if(leptons[0]->isMu() && leptons[1]->isMu()) weight *= 0.706074;
    else weight = 0;
  }
  if(chNum == 129135 || chNum == 147667){
    if(leptons[0]->isEle() && leptons[1]->isEle()) weight *= 0.706074;
    else weight = 0;
  }

  // Trigger
  float trigW = 1;
  if(!m_useMCTrig && includeTrig){
    trigW  = nl == 2 ? m_trigObj->getTriggerWeight(leptons, nt.evt()->isMC, NtSys_NOM) : 1.;
    if(trigW != trigW){ cout<<"\tTrigger weight: "<<trigW<<endl; trigW =0; }// deal with NaN
    if(trigW < 0) trigW = 0;
  }

  // Lepton eff
  /*
  float effW   =  1;
  if(nl > 0 && leptons[0]->isEle){
    Electron* elec = (Electron*) leptons[0];
    effW *= m_susyObj->GetSignalElecSF(elec->clusEta,elec->Pt(),6);
  }
  else if(nl > 0) effW *= leptons[0]->effSF;

  if(nl > 1 && leptons[1]->isEle){
    Electron* elec = (Electron*) leptons[1];
    effW *= m_susyObj->GetSignalElecSF(elec->clusEta,elec->Pt(),6);
  }
  else if(nl > 1) effW *= leptons[1]->effSF;
  */

  float effW   = nl > 0 ? leptons[0]->effSF : 1.;
  effW        *=  nl > 1 ? leptons[1]->effSF : 1.;


  // btag, if included
  float bTag   =  includeBTag ? getBTagWeight(nt.evt()) : 1.;
  return weight * trigW * effW * bTag;
}

/*--------------------------------------------------------------------------------*/
// Get Btag weight
/*--------------------------------------------------------------------------------*/
float SusySelection::getBTagWeight(const Event* evt)
{
  JetVector tempJets;
  for(uint ij=0; ij<m_baseJets.size(); ++ij){
    Jet* jet = m_baseJets.at(ij);
    if( !(jet->Pt() > 20 && fabs(jet->Eta()) < 2.5) ) continue;
    //if( fabs(jet->Eta()) > JET_ETA_CUT  ) continue;
    tempJets.push_back(jet);
  }
  return bTagSF(evt, tempJets, true);
}

/*--------------------------------------------------------------------------------*/
// Photon methods
/*--------------------------------------------------------------------------------*/
bool SusySelection::passCheckMC(int mcRunNumber, float pt)
{

  if( !(159100 <= mcRunNumber && mcRunNumber <= 159104) ) return true;

  // Need to bin the samples in order to combine them
  if( mcRunNumber == 159100 && 25 < pt && pt < 65 ) // Unbinned 17
    return true;
  if( mcRunNumber == 159101 && 65 < pt && pt < 110 ) // Unbinned 35
    return true;
  if( mcRunNumber == 159102 && 110 < pt && pt < 210 ) // Unbinned 70
    return true;
  if( mcRunNumber == 159103 && 210 < pt && pt < 400 ) // Unbinned 140
    return true;
  if(  mcRunNumber == 159104 && 400 < pt             ) // Unbinned 280
    return true;

  // Doesn't match!
  return false;

}
/*--------------------------------------------------------------------------------*/
// Photon xs
/*--------------------------------------------------------------------------------*/
float SusySelection::getPhotonXS(int mcRunNumber)
{

  if( mcRunNumber == 159101)
    return 18652 * 0.57334;
  if( mcRunNumber == 159102)
    return 1640.1 * 0.6398;
  if( mcRunNumber == 159103)
    return 91.747 * 0.77879;
  if( mcRunNumber == 159104)
    return 3.8119 * 0.84244;

  return 1.0;

}

/*--------------------------------------------------------------------------------*/
// Event counters
/*--------------------------------------------------------------------------------*/
void SusySelection::dumpEventCounters()
{


  string v_ET[] = {"ee","mm","em","me"};
  string v_WT[] = {"Raw","Event","Pileup","Pileup A-B3",
		   "LeptonSF","btagSF","TrigSF","All A-B3", "All A-E"};

  for(int w=0; w<WT_N; ++w){
    cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;
    cout << "SusySelection Event counts for weight: " << v_WT[w] << endl;
    cout << endl;
    cout << "read in:       " << n_readin[w]           << endl;
    cout << "pass LAr:      " << n_pass_LAr[w]         << endl;
    cout << "pass BadJet:   " << n_pass_BadJet[w]      << endl;
    cout << "pass BadMu:    " << n_pass_BadMuon[w]     << endl;
    cout << "pass Cosmic:   " << n_pass_Cosmic[w]      << endl;
    cout << "   ------  Start Comparison Here ------ " << endl;
    cout << "pass atleast 2 " << n_pass_atleast2Lep[w] << endl;
    cout << "pass exactly 2 " << n_pass_exactly2Lep[w] << endl;
    cout << "pass nSigLep:  " << n_pass_signalLep[w]   << endl;

    for(int i=0; i<ET_N-1; ++i){
      cout << "************************************" << endl;
      cout << "For dilepton type: " << v_ET[i]       << endl;

      cout << "pass flavor:     " << n_pass_flavor[i][w]    << endl;
      cout << "pass evt trig:   " << n_pass_evtTrig[i][w]   << endl;
      cout << "pass trig match: " << n_pass_trigMatch[i][w] << endl;
      cout << "pass OS:         " << n_pass_os[i][w]        << endl;
      cout << "pass SS:         " << n_pass_ss[i][w]        << endl;
      cout << "-----------------------------------------------------"      << endl;
      cout << "pass SR6 sign:                  " << n_pass_SR6sign[i][w]   << endl;
      cout << "pass SR6 flavor:                " << n_pass_SR6flav[i][w]   << endl;
      cout << "pass SR6 >=2j (no fw veto):     " << n_pass_SR6ge2jNfv[i][w]<< endl;
      cout << "pass SR6 ==2j (no fw veto):     " << n_pass_SR6eq2jNfv[i][w]<< endl;
      cout << "pass SR6 >=1j:                  " << n_pass_SR6ge1j[i][w]   << endl;
      cout << "pass SR6 >=2j:                  " << n_pass_SR6ge2j[i][w]   << endl;
      cout << "pass SR6 ==2j:                  " << n_pass_SR6eq2j[i][w]   << endl;
      cout << "pass SR6 METRel > 50:           " << n_pass_SR6metr[i][w]   << endl;
      cout << "-----------------------------------------------------"      << endl;
      cout << "pass SR7 sign:                  " << n_pass_SR7sign[i][w]   << endl;
      cout << "pass SR7 flavor:                " << n_pass_SR7flav[i][w]   << endl;
      cout << "pass SR7 >=2j (no fw veto):     " << n_pass_SR7ge2jNfv[i][w]<< endl;
      cout << "pass SR7 ==2j (no fw veto):     " << n_pass_SR7eq2jNfv[i][w]<< endl;
      cout << "pass SR7 >=1j:                  " << n_pass_SR7ge1j[i][w]   << endl;
      cout << "pass SR7 >=2j:                  " << n_pass_SR7ge2j[i][w]   << endl;
      cout << "pass SR7 ==2j:                  " << n_pass_SR7eq2j[i][w]   << endl;
      cout << "pass SR7 METRel > 50:           " << n_pass_SR7metr[i][w]   << endl;
      cout << "-----------------------------------------------------"      << endl;
      cout << "pass SR8 sign:                  " << n_pass_SR8sign[i][w]   << endl;
      cout << "pass SR8 flavor:                " << n_pass_SR8flav[i][w]   << endl;
      cout << "pass SR8 >=2j (no fw veto):     " << n_pass_SR8ge2jNfv[i][w]<< endl;
      cout << "pass SR8 ==2j (no fw veto):     " << n_pass_SR8eq2jNfv[i][w]<< endl;
      cout << "pass SR8 >=1j:                  " << n_pass_SR8ge1j[i][w]   << endl;
      cout << "pass SR8 >=2j:                  " << n_pass_SR8ge2j[i][w]   << endl;
      cout << "pass SR8 ==2j:                  " << n_pass_SR8eq2j[i][w]   << endl;
      cout << "pass SR8 METRel > 50:           " << n_pass_SR8metr[i][w]   << endl;
      cout << "-----------------------------------------------------"      << endl;
      cout << "pass SR9 sign:                  " << n_pass_SR9sign[i][w]   << endl;
      cout << "pass SR9 flavor:                " << n_pass_SR9flav[i][w]   << endl;
      cout << "pass SR9 >=2j (no fw veto):     " << n_pass_SR9ge2jNfv[i][w]<< endl;
      cout << "pass SR9 ==2j (no fw veto):     " << n_pass_SR9eq2jNfv[i][w]<< endl;
      cout << "pass SR9 >=1j:                  " << n_pass_SR9ge1j[i][w]   << endl;
      cout << "pass SR9 >=2j:                  " << n_pass_SR9ge2j[i][w]   << endl;
      cout << "pass SR9 ==2j:                  " << n_pass_SR9eq2j[i][w]   << endl;
      cout << "pass SR9 METRel > 50:           " << n_pass_SR9metr[i][w]   << endl;

    }// end loop over event type
  }// end loop over weight type

}

/*--------------------------------------------------------------------------------*/
// Dump Pre objects
/*--------------------------------------------------------------------------------*/
void SusySelection::dumpPreObjects()
{
  ElectronVector preElectrons = getPreElectrons(&nt, NtSys_NOM);
  MuonVector preMuons = getPreMuons(&nt, NtSys_NOM);
  JetVector preJets = getPreJets(&nt, NtSys_NOM);

  cout << "Pre Electrons: "<<preElectrons.size()<<endl;
  for(uint ie=0; ie<preElectrons.size(); ++ie){
    preElectrons[ie]->print();
    //cout<<"Eta: "<<preElectrons[ie]->Eta()<<" Phi: "<<preElectrons[ie]->Phi()<<endl;
  }

  cout<<endl;
  cout << "Pre Muons: "<<preMuons.size() << endl;
  for(uint im=0; im<preMuons.size(); ++im){
    preMuons[im]->print();
    //cout<<"Eta: "<<preMuons[im]->Eta()<<" Phi: "<<preMuons[im]->Phi()<<endl;
  }

  cout << endl;
  cout << "Pre Jets: "<<preJets.size() << endl;
  for(uint ij=0; ij<preJets.size(); ++ij){
    preJets[ij]->print();
    //cout<<"Eta: "<<preJets[ij]->Eta()<<" Phi: "<<preJets[ij]->Phi()<<endl;
  }

}
/*--------------------------------------------------------------------------------*/
// Dump Jets with more information
/*--------------------------------------------------------------------------------*/
void SusySelection::dumpJets()
{

  cout<<"----------------------------------------"<<endl;
  cout<<"Jets:"<<endl;
  for(uint j=0; j<m_signalJets2Lep.size(); ++j){
    const Jet* jet = m_signalJets2Lep.at(j);
    jet->print();
    cout<<"mv1: "<<jet->mv1<<endl;
  }

}

/*--------------------------------------------------------------------------------*/
// Method for checking sysetmatics
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
    for(int j=0; j<2; ++j){
      cout<<"Recalculate phi: "<<j<<endl;
    for(int i=NtSys_SCALEST_UP; i<NtSys_TRIGSF_EL_UP; ++i){
      cout<<"\t"<<i<<endl;
      const Met* tmp_met = getMet(&nt, (SusyNtSys) i,j);
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

  // Now dump the interesting variables
  //--DG-- const Lepton* l0 = leptons[0];
  //--DG-- const Lepton* l1 = leptons[1];
  //--DG-- const Jet*    j  = jets.size() > 0 ? jets[0] : NULL;
  //--DG-- float metRel = getMetRel(met,leptons,jets);

  /*
  out<<"Run "<<nt.evt()->run
     <<" Event "<<nt.evt()->event
     <<" mll "<<Mll(l0,l1)
     <<" l0pt "<<l0->Pt()
     <<" l1pt "<<l1->Pt()
     <<" l0q "<<l0->q
     <<" l1q "<<l1->q
     <<" njets "<<jets.size()
     <<" met "<<met->Et
     <<endl;
  */

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
  for(uint i=0; i<elecs.size(); ++i){
    out<<"Pre Electron: "<<i<<endl;
    printLep((Lepton*) elecs[i]);
  }
  for(uint i=0; i<muons.size(); ++i){
    out<<"Pre Muons: "<<i<<endl;
    printLep((Lepton*) muons[i]);
  }
  for(uint i=0; i<taus.size(); ++i){
    out<<"Pre Taus: "<<i<<endl;
    printLep((Lepton*) taus[i]);
  }
  out<<"++++++++++++++++++++++"<<endl;
  out<<"N Signal Jets"<<jets.size()<<endl;
  for(uint i=0; i<jets.size(); ++i){
    out<<"Jet: "<<i<<endl;
    printJet(jets[i]);
  }
  out<<"++++++++++++++++++++++"<<endl;
  JetVector prejets = getPreJets(&nt, NtSys_NOM);
  out<<"N Pre Jets"<<prejets.size()<<endl;
  for(uint i=0; i<prejets.size(); ++i){
    out<<"Jet: "<<i<<endl;
    printJet(prejets[i]);
  }


  /*
  out<<"Leading Lep:    Pt "<<l0->Pt()<<" Eta "<<l0->Eta()<<" Phi "<<l0->Phi()<<" isComb: "
     <<((Muon*) l0)->isCombined<<" d0 "<<l0->d0<<" d0err "<<l0->errD0
     <<" z0 "<<l0->z0<<" z0err "<<l0->errZ0<<endl;
  out<<"SubLeading Lep: Pt "<<l1->Pt()<<" Eta "<<l1->Eta()<<" Phi "<<l1->Phi()<<" isComb: "
     <<((Muon*) l1)->isCombined<<" d0 "<<l1->d0<<" d0err "<<l1->errD0
     <<" z0 "<<l1->z0<<" z0err "<<l1->errZ0<<endl;
  out<<"Trigger region: "<<m_trigObj->getMMTrigRegion(l0->Pt(),l1->Pt())<<endl;
  out<<"Njets: "<<jets.size()<<endl;
  if(j) out<<"Jet:            Pt "<<j->Pt()<<" Eta "<<j->Eta()<<" Phi "<<j->Phi()<<" jvf: "<<j->jvf<<endl;
  out<<"Met "<<met->Et<<" Phi "<<met->phi<<" metRel "<<metRel<<endl;
  out<<"Mll "<<mll<<" mllj "<<(*l0 + *l1 + *j).M()<<endl;


  // Taus?
  TauVector preTaus = getPreTaus(&nt, NtSys_NOM);
  out<<"Ntaus: "<<preTaus.size()<<endl;
  //out<<"Ntaus: "<<m_baseTaus.size()<<endl;
  //for(uint i=0; i<m_baseTaus.size(); ++i){
  for(uint i=0; i<preTaus.size(); ++i){
    const Tau* tau = m_baseTaus[i];
    out<<"Tau: Pt "<<tau->Pt()<<" Eta "<<tau->Eta()<<" Phi "<< tau->Phi()<<endl;
    out<<"jet BDT "<<tau->jetBDT<<" loose "<<tau->jetBDTSigLoose<<" med "
       <<tau->jetBDTSigMedium<<" tight "<<tau->jetBDTSigTight<<endl;
    out<<"ele BDT "<<tau->eleBDT<<" loose "<<tau->eleBDTLoose<<" med "
       <<tau->eleBDTMedium<<" tight "<<tau->eleBDTTight<<endl;
    out<<"N track: "<<tau->nTrack<<" author "<<tau->author<<endl;
  }

  out<<"Angular inforomation: "<<endl;
  out<<"\tdPhi(l0,l1)   "<<l0->DeltaPhi(*l1)<<endl;
  out<<"\tdPhi(l0,met)  "<<l0->DeltaPhi(met->lv())<<endl;
  out<<"\tdPhi(l1,met)  "<<l1->DeltaPhi(met->lv())<<endl;
  if(j){
    out<<"\tdPhi(l0,jet)  "<<l0->DeltaPhi(*j)<<endl;
    out<<"\tdPhi(l1,jet)  "<<l1->DeltaPhi(*j)<<endl;
    out<<"\tdPhi(met,jet) "<<met->lv().DeltaPhi(*j)<<endl;
    out<<"\tdPhi(ll,jet)  "<<j->DeltaPhi(*l0+*l1)<<endl;
  }
  out<<"\tdPhi(met,ll)  "<<met->lv().DeltaPhi(*l0+*l1)<<endl;
  out<<"MET information: "<<endl;
  out<<"\tRefEle:     "<<met->refEle<<endl;
  out<<"\tRefMuo:     "<<met->refMuo<<endl;
  out<<"\tRefJet:     "<<met->refJet<<endl;
  out<<"\tRefsoftJet: "<<met->softJet<<endl;
  out<<"\tRefGamma:   "<<met->refGamma<<endl;
  out<<"\tRefCell:    "<<met->refCell<<endl;
  out<<endl;
  out<<"Trigger flags: "<<endl;
  out<<"Leading lepton: "<<endl;
  dumpTrigFlag(l0->trigFlags);
  out<<"Subleading lepton: "<<endl;
  dumpTrigFlag(l1->trigFlags);
  */
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
  if(m_dbg) cout << "SusySelection::getXsFromReader()" << endl;
  if(!m_useXsReader || !m_xsReader) return -1.0;
  if(m_xsFromReader < 0.0) {
    m_xsFromReader = m_xsReader->GetXS(  static_cast<int>(nt.evt()->mcChannel) );
  } // end if(xs<0)
  if(m_dbg) cout << "   got " << m_xsFromReader << " for " << static_cast<int>(nt.evt()->mcChannel) << endl;
  return m_xsFromReader;
}
/*--------------------------------------------------------------------------------*/
float SusySelection::computeEventWeightXsFromReader(float lumi)
{
  float defaultXsec = nt.evt()->xsec;
  assert(defaultXsec != 0.0);
  return (getEventWeight(lumi) * getXsFromReader() / defaultXsec);
}
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
