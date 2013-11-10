#include "SusyTest0/SusySelectionMatt.h"

#include <iomanip>
#include "TCanvas.h"

#include "SusyTest0/criteria.h"

using namespace std;
using namespace Susy;

/*--------------------------------------------------------------------------------*/
// SusySelectionMatt Constructor
/*--------------------------------------------------------------------------------*/
SusySelectionMatt::SusySelectionMatt() :
  m_susyObj(NULL),
  m_trigObj(NULL),
  m_useMCTrig(false),
  m_w(1.0),
  m_do1fb(false),
  m_doAD(false),
  m_dumpCounts(true),
  m_nLepMin(2),
  m_nLepMax(2),
  m_cutNBaseLep(true),
  m_ET(ET_Unknown),
  m_doSusy(false),
  m_susyXS(NULL)
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
    n_pass_mll20[w]       = 0;
    n_pass_signalLep[w]   = 0;
    n_pass_HFOR[w]        = 0;
    n_pass_HotSpot[w]     = 0;
    n_pass_TileError[w]   = 0;
    n_pass_FEBCut[w]      = 0;

    // The rest are channel specific.
    for(int i=0; i<ET_N; ++i){
      n_pass_flavor[i][w]   = 0;
      n_pass_signalTau[i][w]   = 0;
      n_pass_evtTrig[i][w]     = 0;
      n_pass_trigMatch[i][w]   = 0;
      n_pass_mll[i][w]         = 0;
      n_pass_ss[i][w]          = 0;
      n_pass_os[i][w]          = 0;
      n_pass_truth[i][w]       = 0;
      // SRmT2a
      n_pass_SRmT2a_zv[i][w] = 0;
      n_pass_SRmT2a_jv[i][w] = 0;
      n_pass_SRmT2a_PtCut[i][w] = 0;
      n_pass_SRmT2a_mt2[i][w] = 0;
      n_pass_SRmT2b_mt2[i][w] = 0;
      n_pass_SRmT2c_mt2[i][w] = 0;
      // SR WW
      n_pass_SRWWa_zv[i][w] = 0;
      n_pass_SRWWa_jv[i][w] = 0;
      n_pass_SRWWa_PtCut[i][w] = 0;
      n_pass_SRWWa_ptll80[i][w] = 0;
      n_pass_SRWWa_metrel80[i][w] = 0;
      n_pass_SRWWa_mll120[i][w] = 0;
      // SR WW b
      n_pass_SRWWb_zv[i][w] = 0;
      n_pass_SRWWb_jv[i][w] = 0;
      n_pass_SRWWb_PtCut[i][w] = 0;
      n_pass_SRWWb_mll170[i][w] = 0;
      n_pass_SRWWb_mt2_90[i][w] = 0;
      // SR WW c
      n_pass_SRWWc_zv[i][w] = 0;
      n_pass_SRWWc_jv[i][w] = 0;
      n_pass_SRWWc_PtCut[i][w] = 0;
      n_pass_SRWWc_mt2_100[i][w] = 0;
      // SR ZJets
      n_pass_SRZjets_zw[i][w] = 0;
      n_pass_SRZjets_2ljets[i][w] = 0;
      n_pass_SRZjets_bfveto[i][w] = 0;
      n_pass_SRZjets_JetPt[i][w] = 0;
      n_pass_SRZjets_PtCut[i][w] = 0;
      n_pass_SRZjets_ptll80[i][w] = 0;
      n_pass_SRZjets_mjjw[i][w] = 0;
      n_pass_SRZjets_metrel80[i][w] = 0;
      n_pass_SRZjets_dRll[i][w] = 0;
      // CR Top mT2
      n_pass_CRTopmt2_OF[i][w] = 0;
      n_pass_CRTopmt2_1bjet[i][w] = 0;
      n_pass_CRTopmt2_lfveto[i][w] = 0;
      n_pass_CRTopmt2_PtCut[i][w] = 0;
      n_pass_CRTopmt2_mt2_70[i][w] = 0;
      // CR Top Met
      n_pass_CRTopMet_OF[i][w] = 0;
      n_pass_CRTopMet_1bjet[i][w] = 0;
      n_pass_CRTopMet_lfveto[i][w] = 0;
      n_pass_CRTopMet_PtCut[i][w] = 0;
      n_pass_CRTopMet_mll120[i][w] = 0;
      n_pass_CRTopMet_ptll80[i][w] = 0;
      n_pass_CRTopMet_metrel80[i][w] = 0;
      // CR WWMet
      n_pass_CRWWMet_OF[i][w] = 0;
      n_pass_CRWWMet_jv[i][w] = 0;
      n_pass_CRWWMet_PtCut[i][w] = 0;
      n_pass_CRWWMet_mll120[i][w] = 0;
      n_pass_CRWWMet_ptll40[i][w] = 0;
      n_pass_CRWWMet_metrel_60_80[i][w] = 0;
      // CR WWmT2
      n_pass_CRWWmt2_OF[i][w] = 0;
      n_pass_CRWWmt2_jv[i][w] = 0;
      n_pass_CRWWmt2_PtCut[i][w] = 0;
      n_pass_CRWWmt2_mt2_50_90[i][w] = 0;
      // CR Top Zjets
      n_pass_CRTopZjets_SF[i][w] = 0;
      n_pass_CRTopZjets_zv[i][w] = 0;
      n_pass_CRTopZjets_2jets[i][w] = 0;
      n_pass_CRTopZjets_bjet[i][w] = 0;
      n_pass_CRTopZjets_fveto[i][w] = 0;
      n_pass_CRTopZjets_PtCut[i][w] = 0;
      n_pass_CRTopZjets_ptll80[i][w] = 0;
      n_pass_CRTopZjets_dRll[i][w] = 0;
      n_pass_CRTopZjets_metrel80[i][w] = 0;
      // CR ZV Met
      n_pass_CRZVMet_zw[i][w] = 0;
      n_pass_CRZVMet_jv[i][w] = 0;
      n_pass_CRZVMet_PtCut[i][w] = 0;
      n_pass_CRZVMet_ptll80[i][w] = 0;
      n_pass_CRZVMet_metrel80[i][w] = 0;
      // CR ZV mt2
      n_pass_CRZVmt2a_zw[i][w] = 0;
      n_pass_CRZVmt2a_jv[i][w] = 0;
      n_pass_CRZVmt2a_PtCut[i][w] = 0;
      n_pass_CRZVmt2a_mt2_90[i][w] = 0;
      n_pass_CRZVmt2b_mt2_120[i][w] = 0;
      n_pass_CRZVmt2c_mt2_150[i][w] = 0;
      n_pass_CRZVmt2d_mt2_100[i][w] = 0;
      // CR ZXZjets
      n_pass_CRZXZjets_SF[i][w] = 0;
      n_pass_CRZXZjets_zw[i][w] = 0;
      n_pass_CRZXZjets_2ljets[i][w] = 0;
      n_pass_CRZXZjets_bfveto[i][w] = 0;
      n_pass_CRZXZjets_JetPtCut[i][w] = 0;
      n_pass_CRZXZjets_PtCut[i][w] = 0;
      n_pass_CRZXZjets_dRll[i][w] = 0;
      n_pass_CRZXZjets_metrel[i][w] = 0;
      n_pass_CRZXZjets_lowerbound[i][w] = 0;

      n_pass_CRWHSS2lss  [i][w] = 0;
      n_pass_CRWHSStauv  [i][w] = 0;
      n_pass_CRWHSSmuiso [i][w] = 0;
      n_pass_CRWHSSeled0 [i][w] = 0;
      n_pass_CRWHSSnfj   [i][w] = 0;
      n_pass_CRWHSSnbj   [i][w] = 0;
      n_pass_CRWHSSnj    [i][w] = 0;
      n_pass_CRWHSS2lpt  [i][w] = 0;
      n_pass_CRWHSSzveto [i][w] = 0;
      n_pass_CRWHSSmwwt  [i][w] = 0;
      n_pass_CRWHSShtmin [i][w] = 0;
      n_pass_CRWHSSmetrel[i][w] = 0;
      n_pass_CRWHSS      [i][w] = 0;
    }
  }// end loop over weight types
  //m_doMuEtconeCut = true;

  m_BoundLow = TF1("expo","expo",1.5,4.0);
  m_BoundLow.FixParameter(0, 6.11861);
  m_BoundLow.FixParameter(1,-1.118815);
}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  if(m_dbg) cout << "SusySelectionMatt::Begin" << endl;

  // Specify 2-lep ana type
  setAnaType(Ana_2LepWH);
  setSelectTaus(true);

  string per = "HCP";
  if(m_do1fb) per = "A-B3";
  if(m_doAD)  per = "A-D7";
  m_trigObj = new DilTrigLogic("Moriond",false/*No Reweight Utils!*/);
  if(m_useMCTrig) m_trigObj->useMCTrigger();

  if(m_doSusy){
    m_susyXS = new XSReader();
    m_susyXS->setDebug(m_dbg);
    m_susyXS->LoadXSInfo();
  }
}

/*--------------------------------------------------------------------------------*/
Bool_t SusySelectionMatt::Process(Long64_t entry)
{
  return kTRUE;
}
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::Terminate()
{
  SusyNtAna::Terminate();
  if(m_dbg) cout << "SusySelectionMatt::Terminate" << endl;
  if(m_dumpCounts)
    dumpEventCounters();
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::selectEvent(bool count)
{
  if(m_dbg) cout << "SusySelectionMatt::selectEvent" << endl;
  // Basic event cuts
  // Upstream Cuts:
  // * GRL
  // * LAr error
  // * tile Error
  // * TTC veto
  // * Primary Veretex
  int flag = nt.evt()->cutFlags[NtSys_NOM];
  // LAr flag (always true..)
  if( !passLAr(flag) )               return false;
  if(count) increment(n_pass_LAr);
  if(m_dbg) cout<<"\tPass LAr"<<endl;
  // Bad Jets
  if( hasBadJet(m_baseJets) )         return false;
  if(count) increment(n_pass_BadJet);
  if(m_dbg) cout<<"\tPass Bad jet"<<endl;
  // FEB issue -- maybe this is calo jets?
  if( !passDeadRegions(m_preJets,
		       m_met,
		       nt.evt()->run,
		       nt.evt()->isMC)
      )                              return false;
  if(count) increment(n_pass_FEBCut);
  if(m_dbg) cout<<"\tPass Dead Regions"<<endl;
  // Bad Muons
  if( hasBadMuon(m_preMuons) )         return false;
  if(count) increment(n_pass_BadMuon);
  if(m_dbg) cout<<"\tPass Bad Muon"<<endl;
  // Cosmic Muons
  if( hasCosmicMuon(m_baseMuons) ) return false;
  if(count) increment(n_pass_Cosmic);
  if(m_dbg) cout<<"\tPass Cosmic"<<endl;
  // HFor
  if( !passHfor() )                  return false;
  if(count) increment(n_pass_HFOR);
  if(m_dbg) cout<<"\tPass Hfor"<<endl;
  // HotSpot
  if( hasHotSpotJet(m_preJets) )       return false;
  if(count) increment(n_pass_HotSpot);
  if(m_dbg) cout<<"\tPass Hot Spot"<<endl;
  // Tile Error
  if( !passTileTripCut(flag) )          return false;
  if(count) increment(n_pass_TileError);
  if(m_dbg) cout<<"\tPass Tile Error"<<endl;
  // If we are not counting (ie doing cutflow) then remove all the baseline muons with eta < 2.5
  if( !count ){
    LeptonVector temp = m_baseLeptons;
    m_baseLeptons.clear();
    for(uint il=0; il<temp.size(); ++il){
      Lepton* lep = temp.at(il);
      if( lep->isMu() && fabs(lep->Eta()) > 2.4) continue;
      m_baseLeptons.push_back(lep);
    }
  }
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::selectBaseEvent(bool doMll, bool count)
{
  // Make sure event level cuts are checked
  if( !selectEvent(count) )               return false;
  // Lepton Analysis cuts
  if( m_baseLeptons.size() < 2 )                return false;
  if(count) increment(n_pass_atleast2Lep);
  if( m_dbg ) cout<<"\tPassed at least 2 leptons"<<endl;
  if(m_baseLeptons.size() != 2)                 return false;
  if(count) increment(n_pass_exactly2Lep);
  if( m_dbg ) cout<<"\tPassed Exctly 2 leptons"<<endl;
  if(doMll && m_baseLeptons.size() == 2){
    if( Mll(m_baseLeptons[0], m_baseLeptons[1]) < 20 )
      return false;
    if(count) increment(n_pass_mll20);
    if( m_dbg ) cout<<"\tPass mll > 20"<<endl;
  }
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::selectAnaEvent(const LeptonVector& leptons, const LeptonVector& baseLeps,
                                       bool count)
{
  if(m_dbg) cout << "SusySelectionMatt::selectAnaEvent" << endl;
  if( !selectBaseEvent(true,count) )           return false;
  m_ET = getDiLepEvtType(baseLeps);
  if(m_ET == ET_me) m_ET = ET_em;
  // Check signal muon if we are counting; otherwise this is handled in the trigger package
  if( count && ( m_ET == ET_em || m_ET == ET_mm) )
    for(uint im=0; im<leptons.size(); ++im)
      if( leptons.at(im)->isMu() && fabs(leptons.at(im)->Eta()) > 2.4 ) return false;
  if( !passNLepCut(leptons) )    return false; // Signal Lepotn Cut
  if( !passTrigger(baseLeps) )   return false; // Trigger Requirement
  if( m_signalTaus.size() != 0 ) return false; // Reject if signal taus
  if(count) increment(n_pass_signalTau[m_ET]);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRmT2a(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // full jet veto
  // pt0 > 35, pt1 > 20
  // metRel > 40
  // Z veto ee/mm
  // mt2 > 90
  if( !oppositeSign(leptons) )          return false;
  if( !passZVeto(leptons) )             return false;
  if(count) increment(n_pass_SRmT2a_zv[m_ET],true,true);
  if( !passJetVeto(jets) )           return false;
  if(count) increment(n_pass_SRmT2a_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )           return false;
  if( leptons[1]->Pt() < 20 )           return false;
  if(count) increment(n_pass_SRmT2a_PtCut[m_ET],true,true);
  if( !loose && getMt2(leptons,met) < 90 )    return false;
  if(count) increment(n_pass_SRmT2a_mt2[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRmT2b(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // full jet veto
  // pt0 > 35, pt1 > 20
  // metRel > 40
  // Z veto ee/mm
  // mt2 > 120
  if( !oppositeSign(leptons) )          return false;
  if( !passZVeto(leptons) )             return false;
  if( !passJetVeto(jets) )           return false;
  if( leptons[0]->Pt() < 35 )           return false;
  if( leptons[1]->Pt() < 20 )           return false;
  if( !loose && getMt2(leptons,met) < 120 )        return false;
  if(count) increment(n_pass_SRmT2b_mt2[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRmT2c(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // full jet veto
  // pt0 > 35, pt1 > 20
  // metRel > 40
  // Z veto ee/mm
  // mt2 > 150
  if( !oppositeSign(leptons) )          return false;
  if( !passZVeto(leptons) )             return false;
  if( !passJetVeto(jets) )           return false;
  if( leptons[0]->Pt() < 35 )           return false;
  if( leptons[1]->Pt() < 20 )           return false;
  if( !loose && getMt2(leptons,met) < 150 )        return false;
  if(count) increment(n_pass_SRmT2c_mt2[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRWWa(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // full jet veto
  // Z veto for ee/mm
  // pt0 > 35, pt1 > 20
  // metRel > 80
  // pt(ll) > 80
  // m(ll) < 120
  if( !oppositeSign(leptons) )              return false;
  if( !passZVeto(leptons) )                    return false;
  if(count) increment(n_pass_SRWWa_zv[m_ET],true,false);
  if( !passJetVeto(jets) )                  return false;
  if(count) increment(n_pass_SRWWa_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_SRWWa_PtCut[m_ET],true,true);
  if( !loose && (*leptons[0]+*leptons[1]).Pt() < 80 ) return false;
  if(count) increment(n_pass_SRWWa_ptll80[m_ET],true,true);
  if( !loose && getMetRel(met,leptons,jets) < 80)     return false;
  if(count) increment(n_pass_SRWWa_metrel80[m_ET],true,true);
  if( Mll(leptons[0],leptons[1]) > 120 )    return false;
  if(count) increment(n_pass_SRWWa_mll120[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRWWb(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // full jet veto
  // Z veto for ee/mm
  // pt0 > 35, pt1 > 20
  // m(ll) < 170
  // mT2 > 100
  if( !oppositeSign(leptons) )              return false;
  if( !passZVeto(leptons) )                    return false;
  if(count) increment(n_pass_SRWWb_zv[m_ET],true);
  if( !passJetVeto(jets) )                  return false;
  if(count) increment(n_pass_SRWWb_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_SRWWb_PtCut[m_ET],true,true);
  if( Mll(leptons[0],leptons[1]) > 170 )    return false;
  if(count) increment(n_pass_SRWWb_mll170[m_ET],true,true);
  if( !loose && getMt2(leptons,met) < 90 )           return false;
  if(count) increment(n_pass_SRWWb_mt2_90[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRWWc(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // full jet veto
  // Z veto
  // pt0 > 35, pt1 > 20
  // mt2 > 100
  if( !oppositeSign(leptons) )              return false;
  if( !passZVeto(leptons) )                    return false;
  if(count) increment(n_pass_SRWWc_zv[m_ET],true);
  if( !passJetVeto(jets) )                  return false;
  if(count) increment(n_pass_SRWWc_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_SRWWc_PtCut[m_ET],true,true);
  if( !loose && getMt2(leptons,met) < 100)         return false;
  if(count) increment(n_pass_SRWWc_mt2_100[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSRZjets(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // SF
  // Z window
  // b and forward veto
  // Require > 1 light jet
  // m(jj) in [50,100] window
  // ptjet1,2 > 45
  // pt(ll) > 80
  // dR(ll) in [0.3,1.5]
  // metRel > 80
  if( !oppositeSign(leptons) )                  return false;
  if( !sameFlavor(leptons) )                    return false;
  if( passZVeto(leptons) )                      return false;
  if(count) increment(n_pass_SRZjets_zw[m_ET],true);
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  if( N_L20 < 2 )                                return false;
  if(count) increment(n_pass_SRZjets_2ljets[m_ET],true,true);
  if( N_B20 + N_F30 != 0 )                       return false;
  if(count) increment(n_pass_SRZjets_bfveto[m_ET],true,true);
  if( jets[0]->Pt() < 45 )                       return false;
  if( jets[1]->Pt() < 45 )                       return false;
  if(count) increment(n_pass_SRZjets_JetPt[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_SRZjets_PtCut[m_ET],true,true);
  if( (*leptons[0]+*leptons[1]).Pt() < 80 )      return false;
  if(count) increment(n_pass_SRZjets_ptll80[m_ET],true,true);
  float mjj = (*jets[0]+*jets[1]).M();
  if( !loose && !(50 < mjj && mjj < 100) )                 return false;
  if(count) increment(n_pass_SRZjets_mjjw[m_ET],true,true);
  if( !loose && getMetRel(met,leptons,jets) < 80)           return false;
  if(count) increment(n_pass_SRZjets_metrel80[m_ET],true,true);
  float dRll = leptons[0]->DeltaR(*leptons[1]);
  if( !loose && !(0.3 < dRll && dRll < 1.5) )              return false;
  if(count) increment(n_pass_SRZjets_dRll[m_ET],true,true);
  return true;
}

/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passPreSRZjets(const LeptonVector& leptons, const JetVector& jets, const Met *met)
{
  // OS
  // SF
  // Z window
  // b and forward veto
  // Require > 1 light jet
  // m(jj) in [50,100] window
  // ptjet1,2 > 45
  // pt(ll) > 80
  // dR(ll) in [0.3,1.5]
  // metRel > 80 <---- Drop metrel
  if( !oppositeSign(leptons) )                  return false;
  if( !sameFlavor(leptons) )                    return false;
  if( passZVeto(leptons) )                      return false;
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  if( N_L20 < 2 )                                return false;
  if( N_B20 + N_F30 != 0 )                       return false;
  if( jets[0]->Pt() < 45 )                       return false;
  if( jets[1]->Pt() < 45 )                       return false;
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if( (*leptons[0]+*leptons[1]).Pt() < 80 )      return false;
  float mjj = (*jets[0]+*jets[1]).M();
  if( !(50 < mjj && mjj < 100) )                 return false;
  float dRll = leptons[0]->DeltaR(*leptons[1]);
  if( !(0.3 < dRll && dRll < 1.5) )              return false;
  if( !nt.evt()->isMC && getMetRel(met,leptons,jets) > 80)  return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZjets(const LeptonVector& leptons, const JetVector& jets, const Met *met)
{
  // Standard
  if( !oppositeSign(leptons) )                  return false;
  if( !sameFlavor(leptons) )                    return false;
  if( passZVeto(leptons) )                      return false;
  // Selection of jets
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  if( N_L20 < 2 )                                return false;
  if( N_B20 + N_F30 != 0 )                       return false;
  // Pt of jets
  if( jets[0]->Pt() < 45 )                       return false;
  if( jets[1]->Pt() < 45 )                       return false;
  // Pt of leptons
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  // Reverse Pt(ll)
  if( (*leptons[0]+*leptons[1]).Pt() > 80 )      return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passVRSS(const LeptonVector& leptons, const JetVector& jets, const Met *met)
{
  if(m_dbg) cout << "SusySelectionMatt::passVRSS" << endl;
  // Updating Validation Region for new SS SR
  // * Same sign
  // * EtRel > 40 GeV
  // * ee: Z veto
  // * mm: 90 < m(ll) < 120
  // Only SS
  if( !sameSign(leptons) )               return false;
  // MetRel
  float metRel = getMetRel(met,leptons,jets);
  if( metRel < 40 )                      return false;
  // m(ll) selection
  DiLepEvtType et = getDiLepEvtType(leptons);
  if( et == ET_ee && !passZVeto(leptons) ) return false;
  if( et == ET_mm ){
    float mll = Mll(leptons[0],leptons[1]);
    if( 90 < mll && mll < 120 )            return false;
  }
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRWWMet(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // jet veto
  // pt0 > 35, pt1 > 20
  // 60 < metRel < 80
  // ptll > 40
  if( !oppositeSign(leptons) )                return false;
  if( !oppositeFlavor(leptons) )              return false;
  if(count) increment(n_pass_CRWWMet_OF[m_ET],true);
  if( !passJetVeto(jets) )                    return false;
  if(count) increment(n_pass_CRWWMet_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )                 return false;
  if( leptons[1]->Pt() < 20 )                 return false;
  if(count) increment(n_pass_CRWWMet_PtCut[m_ET],true,true);
  if( Mll(leptons[0],leptons[1]) >120 )       return false;
  if(count) increment(n_pass_CRWWMet_mll120[m_ET],true,true);
  if( (*leptons[0]+*leptons[1]).Pt() < 40 )   return false;
  if(count) increment(n_pass_CRWWMet_ptll40[m_ET],true,true);
  float metRel = getMetRel(met,leptons,jets);
  if( !loose && !(60 < metRel && metRel < 80) )         return false;
  if(count) increment(n_pass_CRWWMet_metrel_60_80[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRWWmT2(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // Jet Veto
  // pt0 > 35, pt1>20
  // 50 < mt2 < 90
  if( !oppositeSign(leptons) )              return false;
  if( !oppositeFlavor(leptons) )            return false;
  if(count) increment(n_pass_CRWWmt2_OF[m_ET],true);
  if( !passJetVeto(jets) )                  return false;
  if(count) increment(n_pass_CRWWmt2_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_CRWWmt2_PtCut[m_ET],true,true);
  float mT2 = getMt2(leptons,met);
  if( !loose && !(50 < mT2 && mT2 < 90) )         return false;
  if(count) increment(n_pass_CRWWmt2_mt2_50_90[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRTopMet(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // B20 > 0
  // pt0 > 35, pt1 > 20
  // metRel > 80
  // pt(ll) > 80
  // m(ll) < 120
  if( !oppositeSign(leptons) )              return false;
  if( !oppositeFlavor(leptons) )            return false;
  if(count) increment(n_pass_CRTopMet_OF[m_ET],true);
  if( numberOfCBJets(jets) == 0 )           return false;
  if(count) increment(n_pass_CRTopMet_1bjet[m_ET],true,true);
  if( numberOfCLJets(jets) != 0 )           return false;
  if( numberOfFJets(jets) != 0 )            return false;
  if(count) increment(n_pass_CRTopMet_lfveto[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_CRTopMet_PtCut[m_ET],true,true);
  if( Mll(leptons[0],leptons[1]) > 120 )    return false;
  if(count) increment(n_pass_CRTopMet_mll120[m_ET],true,true);
  if( (*leptons[0]+*leptons[1]).Pt() < 80 ) return false;
  if(count) increment(n_pass_CRTopMet_ptll80[m_ET],true,true);
  if( !loose && getMetRel(met,leptons,jets) < 80 )    return false;
  if(count) increment(n_pass_CRTopMet_metrel80[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRTopmT2(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // Top Control Region
  if( !oppositeSign(leptons) )              return false;
  if( !oppositeFlavor(leptons) )            return false;
  if(count) increment(n_pass_CRTopmt2_OF[m_ET],true,false);
  int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  if( N_B20 == 0 )                           return false;
  if(count) increment(n_pass_CRTopmt2_1bjet[m_ET],true,true);
  if( N_L20 + N_F30 !=0 )                   return false;
  if(count) increment(n_pass_CRTopmt2_lfveto[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )                  return false;
  if( leptons[1]->Pt() < 20 )                  return false;
  if(count) increment(n_pass_CRTopmt2_PtCut[m_ET],true,true);
  if( !loose && getMt2(leptons,met) < 70 )               return false;
  if(count) increment(n_pass_CRTopmt2_mt2_70[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRTopZjets(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // SF
  // Z Veto
  // metRel > 80
  // L20 + B20 >= 2
  // B20 >=1
  // F30 == 0
  // Asymmetric Pt Cut
  // pt(ll) > 80
  // 0.3 < dR(ll) < 1.5
  if( !oppositeSign(leptons) )              return false;
  if( !sameFlavor(leptons) )                return false;
  if(count) increment(n_pass_CRTopZjets_SF[m_ET],true);
  if( !passZVeto(leptons) )                 return false;
  if(count) increment(n_pass_CRTopZjets_zv[m_ET],true);
  int N_B20 = numberOfCBJets(jets);
  if(count) increment(n_pass_CRTopZjets_2jets[m_ET],true,true);
  if( N_B20 == 0 )                          return false;
  if(count) increment(n_pass_CRTopZjets_bjet[m_ET],true,true);
  if( numberOfFJets(jets) != 0 )            return false;
  if(count) increment(n_pass_CRTopZjets_fveto[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )               return false;
  if( leptons[1]->Pt() < 20 )               return false;
  if(count) increment(n_pass_CRTopZjets_PtCut[m_ET],true,true);
  if( (*leptons[0]+*leptons[1]).Pt() < 80 ) return false;
  if(count) increment(n_pass_CRTopZjets_ptll80[m_ET],true,true);
  float dRll = leptons[0]->DeltaR(*leptons[1]);
  if( !(0.3 < dRll && dRll < 1.5) )         return false;
  if(count) increment(n_pass_CRTopZjets_dRll[m_ET],true,true);
  if( !loose && getMetRel(met,leptons,jets) < 80 )    return false;
  if(count) increment(n_pass_CRTopZjets_metrel80[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZXZjets(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS SF
  // Z selection
  // Asymmetric pt cut
  // L20 >= 2
  // Veto b and forward jets
  // dR(ll) in [1.4,4.0]
  // MetRel > 80 GeV
  // BoundLow(dR(ll)) <--- requires a fit funciton
  if( !oppositeSign(leptons) )               return false;
  if( !sameFlavor(leptons) )                 return false;
  if(count) increment(n_pass_CRZXZjets_SF[m_ET],true);
  if( passZVeto(leptons) )                   return false;
  if(count) increment(n_pass_CRZXZjets_zw[m_ET],true);
  if( numberOfCLJets(jets) < 2 )             return false;
  if(count) increment(n_pass_CRZXZjets_2ljets[m_ET],true,true);
  if( numberOfCBJets(jets) != 0 )             return false;
  if( numberOfFJets(jets) != 0 )              return false;
  if(count) increment(n_pass_CRZXZjets_bfveto[m_ET],true,true);
  if( jets[0]->Pt() < 35 )                    return false;
  if( jets[1]->Pt() < 25 )                    return false;
  if(count) increment(n_pass_CRZXZjets_JetPtCut[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )                 return false;
  if( leptons[1]->Pt() < 20 )                 return false;
  if(count) increment(n_pass_CRZXZjets_PtCut[m_ET],true,true);
  float dRll = leptons[0]->DeltaR(*leptons[1]);
  if( !(1.5 < dRll && dRll < 4.0) )           return false;
  if(count) increment(n_pass_CRZXZjets_dRll[m_ET],true,true);
  if( !loose && getMetRel(met,leptons,jets) < 80 )      return false;
  if(count) increment(n_pass_CRZXZjets_metrel[m_ET],true,true);
  float ptll = (*leptons[0]+*leptons[1]).Pt();
  if( ptll < m_BoundLow.Eval(dRll) )         return false;
  if(count) increment(n_pass_CRZXZjets_lowerbound[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZVMet(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // Jet Veto
  // Z window
  // pt0 > 35, pt1 > 20
  // metRel > 80
  // pt(ll) > 80
  if( !oppositeSign(leptons) )               return false;
  if( passZVeto(leptons) )                   return false;
  if(count) increment(n_pass_CRZVMet_zw[m_ET],true);
  if( !passJetVeto(jets) )                   return false;
  if(count) increment(n_pass_CRZVMet_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )                return false;
  if( leptons[1]->Pt() < 20 )                return false;
  if(count) increment(n_pass_CRZVMet_PtCut[m_ET],true,true);
  if( !loose && (*leptons[0]+*leptons[1]).Pt()  < 80 ) return false;
  if(count) increment(n_pass_CRZVMet_ptll80[m_ET],true,true);
  if( !loose && getMetRel(met,leptons,jets) < 80 )     return false;
  if(count) increment(n_pass_CRZVMet_metrel80[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZVmT2a(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // jet veto
  // Z window
  // pt0 > 35, pt1 > 20
  // mt2 > 90
  if( !oppositeSign(leptons) )   return false;
  if( passZVeto(leptons) )       return false;
  if(count) increment(n_pass_CRZVmt2a_zw[m_ET],true);
  if( !passJetVeto(jets) )       return false;
  if(count) increment(n_pass_CRZVmt2a_jv[m_ET],true,true);
  if( leptons[0]->Pt() < 35 )    return false;
  if( leptons[1]->Pt() < 20 )    return false;
  if(count) increment(n_pass_CRZVmt2a_PtCut[m_ET],true,true);
  if( !loose && getMt2(leptons,met) < 90 ) return false;
  if(count) increment(n_pass_CRZVmt2a_mt2_90[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZVmT2b(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // jet veto
  // Z window
  // pt0 > 35, pt1 > 20
  // mt2 > 90
  if( !oppositeSign(leptons) )   return false;
  if( passZVeto(leptons) )       return false;
  if( !passJetVeto(jets) )       return false;
  if( leptons[0]->Pt() < 35 )    return false;
  if( leptons[1]->Pt() < 20 )    return false;
  if( !loose && getMt2(leptons,met) < 120 ) return false;
  if(count) increment(n_pass_CRZVmt2b_mt2_120[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZVmT2c(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // jet veto
  // Z window
  // pt0 > 35, pt1 > 20
  // mt2 > 90
  if( !oppositeSign(leptons) )   return false;
  if( passZVeto(leptons) )       return false;
  if( !passJetVeto(jets) )       return false;
  if( leptons[0]->Pt() < 35 )    return false;
  if( leptons[1]->Pt() < 20 )    return false;
  if( !loose && getMt2(leptons,met) < 150 ) return false;
  if(count) increment(n_pass_CRZVmt2c_mt2_150[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRZVmT2d(const LeptonVector& leptons, const JetVector& jets, const Met *met, bool count, bool loose)
{
  // OS
  // jet veto
  // Z window
  // pt0 > 35, pt1 > 20
  // mt2 > 90
  if( !oppositeSign(leptons) )   return false;
  if( passZVeto(leptons) )       return false;
  if( !passJetVeto(jets) )       return false;
  if( leptons[0]->Pt() < 35 )    return false;
  if( leptons[1]->Pt() < 20 )    return false;
  if( !loose && getMt2(leptons,met) < 100 ) return false;
  if(count) increment(n_pass_CRZVmt2d_mt2_100[m_ET],true,true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCRPremT2(const LeptonVector& leptons, const JetVector& jets, const Met *met)
{
  // OS
  // full jet veto
  // pt0 > 35, pt1 > 20
  // metRel > 40
  // Z veto ee/mm
  if( !oppositeSign(leptons) )          return false;
  if( !passZVeto(leptons) )             return false;
  if( !passJetVeto(jets) )           return false;
  if( leptons[0]->Pt() < 35 )           return false;
  if( leptons[1]->Pt() < 20 )           return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSimpleZ(const LeptonVector& leptons, const JetVector& jets, const Met *met)
{
  // OS
  // Z window
  // Met > 40
  if( sameSign(leptons) ) return false;
  if( passZVeto(leptons) ) return false;
  float metrel = getMetRel(met, leptons, jets);
  if( metrel < 40 ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passSimpleZ2(const LeptonVector& leptons, const JetVector& jets, const Met *met)
{
  // OS
  // Z window
  // 70 < met < 100
  // jet veto
  if( sameSign(leptons) ) return false;
  if( passZVeto(leptons) ) return false;
  float metrel = getMetRel(met,leptons,jets);
  if( !(70 < metrel && metrel < 100) ) return false;
  if( !passJetVeto(jets) ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passMuonRelIso(const LeptonVector &leptons, float maxVal)
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
/*--------------------------------------------------------------------------------*/
SsPassFlags SusySelectionMatt::passWhSS(const LeptonVector& leptons, const JetVector& jets, const Met* met)
{
  // for now keep it simple:
  // - 2lep ss (I assume we don't need qFlip here)
  // - pt0>30
  // -  pt1>20 for em, ee
  // - ht>200
  // - d0 < 3 for el
  // - z veto (10GeV) for ee
  // - mww > 150 (ee), > 140 (em), >100 (mm)
  // - metrel > 50 for ee, em
  //
  // DG: be careful, this is slightly different from SusySelection::passSrSs.
  // In particular, here we don't require the trigger match and some
  // of the criteria are in a different order (e.g. same-sign).
  // At some point try to unify SusySelection with SusySelectionMatt.
  SsPassFlags f;
  bool lsf(false), bsf(false); // compute trigw and btagw only when accepting the event
  if(sameSign(leptons)) { increment(n_pass_CRWHSS2lss  [m_ET], lsf, bsf); f.sameSign=true;} else return f;
  DiLepEvtType ll = m_ET = getDiLepEvtType(leptons);
  bool isee(ll==ET_ee), isem(ll==ET_em||ll==ET_me), ismm(ll==ET_mm);
  float ptL0Min  = 30;
  float ptL1Min  = (ismm ? 0.0 : 20.0);
  float htMin    = 200;
  bool applyMllZveto(isee);
  float mZ0(91.2);
  float loMllZ(applyMllZveto ? mZ0-10. : FLT_MAX);
  float hiMllZ(applyMllZveto ? mZ0+10. : FLT_MIN);
  float mtwwMin = (isee ? 150 : (isem ? 140 : (ismm ? 100 : FLT_MIN))); // todo : for now keep it simple, just one cut
  float metRelMin = (isee ? 50 : (isem ? 50 : (ismm ? FLT_MIN : FLT_MIN))); // for now simple

  if(m_signalTaus.size()==0)                          { increment(n_pass_CRWHSStauv  [m_ET], lsf, bsf); f.tauVeto=true;} else  return f;
  if(numberOfFJets(jets)==0)                          { increment(n_pass_CRWHSSnfj   [m_ET], lsf, bsf); f.fjveto =true;} else  return f;
  if(numberOfCBJets(jets)==0)                         { increment(n_pass_CRWHSSnbj   [m_ET], lsf, bsf); f.bjveto =true;} else  return f;
  if(numberOfCLJets(jets)>0)                          { increment(n_pass_CRWHSSnj    [m_ET], lsf, bsf); f.ge1j   =true;} else  return f;
  if(susy::pass2LepPt    (leptons, ptL0Min, ptL1Min)) { increment(n_pass_CRWHSS2lpt  [m_ET], lsf, bsf); f.lepPt  =true;} else  return f;
  if(susy::passZllVeto   (leptons, loMllZ, hiMllZ))   { increment(n_pass_CRWHSSzveto [m_ET], lsf, bsf); f.zllVeto=true;} else  return f;
  if(susy::passMtLlMetMin(leptons, met, mtwwMin))     { increment(n_pass_CRWHSSmwwt  [m_ET], lsf, bsf); f.mtllmet=true;} else  return f;
  if(susy::passHtMin     (leptons, jets, met, htMin)) { increment(n_pass_CRWHSShtmin [m_ET], lsf, bsf); f.ht     =true;} else  return f;
  if(getMetRel(met,leptons,jets)>metRelMin)           { increment(n_pass_CRWHSSmetrel[m_ET], lsf, bsf); f.metrel =true;} else  return f;
  lsf = bsf = true;
  increment(n_pass_CRWHSS[m_ET], lsf, bsf);
  return f;
}

/*--------------------------------------------------------------------------------*/
// Generic cuts
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passHfor()
{
  if(nt.evt()->hfor == 4 ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passNLepCut(const LeptonVector& leptons)
{
  uint nLep = leptons.size();
  if(m_nLepMin>=0 && nLep < m_nLepMin) return false;
  if(m_nLepMax>=0 && nLep > m_nLepMax) return false;
  increment(n_pass_signalLep); //+=m_w;
  increment(n_pass_flavor[m_ET],true);
  // To stay inline with Anders' cutflow
  if(m_ET == ET_me) increment(n_pass_flavor[ET_em],true);
  if(m_ET == ET_em) increment(n_pass_flavor[ET_me],true);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passNBaseLepCut(const LeptonVector& baseLeptons)
{
  if(m_cutNBaseLep){
    uint nLep = baseLeptons.size();
    if(m_nLepMin>=0 && nLep < m_nLepMin) return false;
    if(m_nLepMax>=0 && nLep > m_nLepMax) return false;
  }
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passTrigger(const LeptonVector& leptons)
{
  if(leptons.size() != 2)
    return false;
  bool passEvtTrig   = m_trigObj->passDilEvtTrig(leptons, m_met->Et, nt.evt());
  bool passTrigMatch = m_trigObj->passDilTrigMatch(leptons, m_met->Et, nt.evt());
  if( passEvtTrig ){
    increment(n_pass_evtTrig[m_ET],true);
  }
  if( passEvtTrig && passTrigMatch){
    increment(n_pass_trigMatch[m_ET],true);
    return true;
  }
  return false;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::sameFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return (leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::oppositeFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return !(leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::sameSign(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return leptons.at(0)->q * leptons.at(1)->q > 0;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::oppositeSign(const LeptonVector& leptons)
{
  return !(sameSign(leptons));
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passMll(const LeptonVector& leptons, float mll)
{
  if( leptons.size() < 2 ) return false;
  if( (*leptons.at(0) + *leptons.at(1)).M() < mll ) return false;
  increment(n_pass_mll[m_ET]);
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passBadMet(const Met* met, float cutval)
{
  JetVector preJets = getPreJets(&nt, NtSys_NOM);
  float RefEle = met->refEle;
  float px     = met->refJet_etx;
  float py     = met->refJet_ety;
  for(uint ij=0; ij<preJets.size(); ++ij){
    px += preJets.at(ij)->Px();
    py += preJets.at(ij)->Py();
  }
  float OJ_RefEle = sqrt(px*px+py*py)/ RefEle;
  return OJ_RefEle > cutval;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passJetVeto(const JetVector& jets)
{
  // Require no light, b, or forward jets
  int N_L25 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L25 + N_B20 + N_F30 == 0);
}
/*--------------------------------------------------------------------------------*/
int SusySelectionMatt::nL20Close(const JetVector& jets)
{
  // number of L20 jets with dR > 1.0 to leading b jet
  Jet* leadingBjet = NULL;
  int nL20 = 0;
  for(uint ij=0; ij<jets.size(); ++ij){
    if( !isCentralBJet(jets.at(ij)) ) continue;
    if( !leadingBjet ) leadingBjet = jets.at(ij);
    else if(leadingBjet->Pt() < jets.at(ij)->Pt())
      leadingBjet = jets.at(ij);
  }// end loop
  if( !leadingBjet ) return nL20;
  for(uint ij=0; ij<jets.size(); ++ij){
    if( !isCentralLightJet(jets.at(ij)) ) continue;
    if( leadingBjet->DeltaR( *jets.at(ij) ) > 1.0 ) nL20++;
  }
  return nL20;
}
/*--------------------------------------------------------------------------------*/
int SusySelectionMatt::nB20Close(const JetVector& jets)
{
  // number of L20 jets with dR > 1.0 to leading b jet
  Jet* leadingBjet = NULL;
  int nB20 = 0;
  for(uint ij=0; ij<jets.size(); ++ij){
    if( !isCentralBJet(jets.at(ij)) ) continue;
    if( !leadingBjet ) leadingBjet = jets.at(ij);
    else if(leadingBjet->Pt() < jets.at(ij)->Pt())
      leadingBjet = jets.at(ij);
  }// end loop
  if( !leadingBjet ) return nB20;
  for(uint ij=0; ij<jets.size(); ++ij){
    if( !isCentralBJet(jets.at(ij)) ) continue;
    if( leadingBjet == jets.at(ij) )      continue;
    if( leadingBjet->DeltaR( *jets.at(ij) ) > 1.0 ) nB20++;
  }
  return nB20;
}
/*--------------------------------------------------------------------------------*/
int SusySelectionMatt::nF30Close(const JetVector& jets)
{
  // number of L20 jets with dR > 1.0 to leading b jet
  Jet* leadingBjet = NULL;
  int nF30 = 0;
  for(uint ij=0; ij<jets.size(); ++ij){
    if( !isCentralBJet(jets.at(ij)) ) continue;
    if( !leadingBjet ) leadingBjet = jets.at(ij);
    else if(leadingBjet->Pt() < jets.at(ij)->Pt())
      leadingBjet = jets.at(ij);
  }// end loop
  if( !leadingBjet ) return nF30;
  for(uint ij=0; ij<jets.size(); ++ij){
    if( !isForwardJet(jets.at(ij)) ) continue;
    if( leadingBjet->DeltaR( *jets.at(ij) ) > 1.0 ) nF30++;
  }
  return nF30;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passbJetVeto(const JetVector& jets)
{
  // Reject if there is a b jet using 2L definition
  int N_B20 = numberOfCBJets(jets);
  return (N_B20 == 0);
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passge2Jet(const JetVector& jets)
{
  int N_L25 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);
  return (N_L25 >=2 && N_B20 + N_F30 == 0);
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passZVeto(const LeptonVector& leptons, float Zlow, float Zhigh)
{
  if( leptons.size() < 2 )   return false;
  if( !sameFlavor(leptons) ) return true;
  float mll = (*leptons.at(0) + *leptons.at(1)).M();
  if( Zlow < mll && mll < Zhigh ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passMETRel(const Met *met, const LeptonVector& leptons,
				 const JetVector& jets, float metMax){
  if( getMetRel(met,leptons,jets) < metMax ) return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passdPhi(TLorentzVector v0, TLorentzVector v1, float cut)
{
  return v0.DeltaPhi(v1) > cut;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passMT2(const LeptonVector& leptons, const Met* met, float cut)
{
  float mt2 = getMt2(leptons, met);
  return (mt2 > cut);
}
/*--------------------------------------------------------------------------------*/
float SusySelectionMatt::getMt2(const LeptonVector& leptons, const Met* met)
{
  if( leptons.size() < 2 ) return false;
  TLorentzVector metlv = met->lv();
  TLorentzVector l0    = *leptons.at(0);
  TLorentzVector l1    = *leptons.at(1);
  double pTMiss[3] = {0.0, metlv.Px(), metlv.Py()};
  double pA[3]     = {0.0, l0.Px(), l0.Py()};
  double pB[3]     = {0.0, l1.Px(), l1.Py()};
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pA,pB,pTMiss);
  mt2_event.set_mn(0); // LSP mass = 0 is Generic
  return mt2_event.get_mt2();
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isRealLepton(const Lepton* lep)
{
  // Updated way of handling real and fake leptons using LeptonTruthTools
  // Need to handle new g2ww -- Assume all real for now
  uint mcId = nt.evt()->mcChannel;
  if( mcId == 169471 || mcId == 169472 || mcId == 169473 || mcId == 169474 ||
      mcId == 169475 || mcId == 169476 || mcId == 169477 || mcId == 169478 ||
      mcId == 169479)
    return true;
  return (lep->truthType == RecoTruthMatch::PROMPT);
  // Code taken from Steve.  There seems to be an issue with Sherpa samples, so
  // need to handle those separately. Also just for clarification:
  // * mcOrigin = 9 -- Tau Lepton
  // * mcType   = 1 -- Unknown Electron
  // * mcType   = 2 -- Iso Electron
  // * mcType   = 5 -- Unknown Muon
  // * mcType   = 6 -- Iso Muon
  // Cut is sample dependent due to Sherpa classifications broken
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
bool SusySelectionMatt::isFakeLepton(const Lepton* lep)
{
  return !isRealLepton(lep);
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isConvLepton(const Lepton* lep)
{
  //return lep->mcOrigin == 5;
  bool isConv       = lep->truthType == RecoTruthMatch::CONV;
  //bool isConv       = lep->mcOrigin == 5;
  bool isChargeFlip =  lep->isEle() ? ((Electron*) lep)->isChargeFlip : false;
  return isConv && !isChargeFlip;
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isHFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::HF);
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isLFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::LF);
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isQCDLepton(const Lepton* lep)
{
  return isHFLepton(lep) || isLFLepton(lep);
}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isTrueDilepton(const LeptonVector &leptons)
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
bool SusySelectionMatt::isFakeDilepton(const LeptonVector &leptons)
{
  // Maybe not 100% kosher, but I just want to make sure I am not
  // double counting, so I require dilepton events to be real
  if( leptons.size() != 2 ) return false;
  bool l0_fake = isFakeLepton(leptons[0]);
  bool l1_fake = isFakeLepton(leptons[1]);
  bool l0_cf = leptons[0]->isEle() ? ((Electron*) leptons[0])->isChargeFlip : false;
  bool l1_cf = leptons[1]->isEle() ? ((Electron*) leptons[1])->isChargeFlip : false;
  return (l0_fake || l1_fake) && !l0_cf && !l1_cf; // ignoring charge flip
}
/*--------------------------------------------------------------------------------*/
float SusySelectionMatt::getEvtWeight(const LeptonVector& leptons, bool includeBTag, bool includeTrig,
				  bool doMediumpp)
{
  if( !nt.evt()->isMC ) return 1.;
  uint nl = leptons.size();
  float weight = 1;
  // lumi, xs, sumw, pileup
  if(m_do1fb) weight = getEventWeightAB3();
  else if(m_doAD)  weight = getEventWeightFixed(nt.evt()->mcChannel, LUMI_A_D);
  else if(m_doSusy){
    Event* evt = nt.evt();
    float weight = evt->w * evt->wPileup * LUMI_A_E / evt->sumw;
    weight *= m_susyXS->GetXS( evt->mcChannel );
  }
  else weight = getEventWeight(LUMI_A_L, true);
  // Trigger
  float trigW = 1;
  if(!m_useMCTrig && includeTrig){
    trigW  = nl == 2 ? m_trigObj->getTriggerWeight(leptons,
						   nt.evt()->isMC,
						   m_met->Et,
						   m_signalJets2Lep.size(),
						   nt.evt()->nVtx,
						   NtSys_NOM) : 1.;
    if(trigW != trigW){ cout<<"\tTrigger weight: "<<trigW<<endl; trigW =0; }// deal with NaN
    if(trigW < 0) trigW = 0;
  }
  float effW   = nl > 0 ? leptons[0]->effSF : 1.;
  effW        *=  nl > 1 ? leptons[1]->effSF : 1.;
  // btag, if included
  float bTag   =  includeBTag ? getBTagWeight(nt.evt()) : 1.;
  return weight * trigW * effW * bTag;
}
/*--------------------------------------------------------------------------------*/
float SusySelectionMatt::getBTagWeight(const Event* evt)
{
  JetVector tempJets;
  for(uint ij=0; ij<m_baseJets.size(); ++ij){
    Jet* jet = m_baseJets.at(ij);
    if( !(jet->Pt() > 20 && fabs(jet->detEta) < 2.4) ) continue;
    tempJets.push_back(jet);
  }
  return bTagSF(evt, tempJets, evt->mcChannel, BTag_NOM);
}
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::dumpEventCounters()
{
  string v_WT[] = {"Raw","Event","Pileup","Pileup A-B3",
                   "LeptonSF","btagSF","TrigSF","All A-B3", "All A-E"};

  for(int w=0; w<WT_N; ++w){
    cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;
    cout << "SusySelectionMatt Event counts for weight: " << v_WT[w] << endl;
    cout << endl;
    cout << "read in:       " << n_readin[w]           << endl;
    cout << "pass LAr:      " << n_pass_LAr[w]         << endl;
    cout << "pass BadJet:   " << n_pass_BadJet[w]      << endl;
    cout << "pass FEB:      " << n_pass_FEBCut[w]      << endl;
    cout << "pass BadMu:    " << n_pass_BadMuon[w]     << endl;
    cout << "pass Cosmic:   " << n_pass_Cosmic[w]      << endl;
    cout << "pass HFOR:     " << n_pass_HFOR[w]        << endl;
    cout << "pass Hot Spot: " << n_pass_HotSpot[w]     << endl;
    cout << "pass Tile Err: " << n_pass_TileError[w]   << endl;
    cout << "   ------  Start Comparison Here ------ " << endl;
    cout << ">=2 Lep        " << n_pass_atleast2Lep[w] << endl;
    cout << "pass Exactly 2 " << n_pass_exactly2Lep[w] << endl;
    cout << "pass mll > 20  " << n_pass_mll20[w]       << endl;
    cout << "pass nSigLep:  " << n_pass_signalLep[w]   << endl;

    cout << "************************************" << endl;

     cout << "Cut                   \tee\t\tmm\t\tem" << endl;
    printCounter("pass flavor:     " , n_pass_flavor    , w);
    printCounter("pass evt trig:   " , n_pass_evtTrig   , w);
    printCounter("pass trig match: " , n_pass_trigMatch , w);
    printCounter("pass tau veto:   " , n_pass_signalTau , w);
    printCounter("pass Truth:      " , n_pass_truth     , w);
    printCounter("pass OS:         " , n_pass_os        , w);
    printCounter("pass SS:         " , n_pass_ss        , w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRmT2 ZVeto            ", n_pass_SRmT2a_zv, w);
    printCounter("pass: SRmT2 Jet Veto         ", n_pass_SRmT2a_jv, w);
    printCounter("pass: SRmT2 PtCut            ", n_pass_SRmT2a_PtCut, w);
    printCounter("pass: SRmT2 mT2 > 90         ", n_pass_SRmT2a_mt2, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRmT2 mT2 > 120        ", n_pass_SRmT2b_mt2, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRmT2 mT2 > 150        ", n_pass_SRmT2c_mt2, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRWWa ZVeto            ", n_pass_SRWWa_zv, w);
    printCounter("pass: SRWWa JetVeto          ", n_pass_SRWWa_jv, w);
    printCounter("pass: SRWWa PtCut            ", n_pass_SRWWa_PtCut, w);
    printCounter("pass: SRWWa pt(ll) > 80      ", n_pass_SRWWa_ptll80, w);
    printCounter("pass: SRWWa metrel > 80      ", n_pass_SRWWa_metrel80, w);
    printCounter("pass: SRWWa m(ll) < 120      ", n_pass_SRWWa_mll120, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRWWb ZVeto            ", n_pass_SRWWb_zv, w);
    printCounter("pass: SRWWb JetVeto          ", n_pass_SRWWb_jv, w);
    printCounter("pass: SRWWb PtCut            ", n_pass_SRWWb_PtCut, w);
    printCounter("pass: SRWWb m(ll) < 170      ", n_pass_SRWWb_mll170, w);
    printCounter("pass: SRWWb mT2 > 90         ", n_pass_SRWWb_mt2_90, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRWWc ZVeto            ", n_pass_SRWWc_zv, w);
    printCounter("pass: SRWWc JetVeto          ", n_pass_SRWWc_jv, w);
    printCounter("pass: SRWWc PtCut            ", n_pass_SRWWc_PtCut, w);
    printCounter("pass: SRWWc mT2 > 100        ", n_pass_SRWWc_mt2_100, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: SRZjets Z window       ", n_pass_SRZjets_zw, w);
    printCounter("pass: SRZjets >= 1 LJets     ", n_pass_SRZjets_2ljets, w);
    printCounter("pass: SRZjets forward/b veto ", n_pass_SRZjets_bfveto, w);
    printCounter("pass: SRZjets Jet0,1 > 45    ", n_pass_SRZjets_JetPt, w);
    printCounter("pass: SRZjets PtCut          ", n_pass_SRZjets_PtCut, w);
    printCounter("pass: SRZjets pt(ll) > 80    ", n_pass_SRZjets_ptll80, w);
    printCounter("pass: SRZjets 50<m(jj)<100   ", n_pass_SRZjets_mjjw, w);
    printCounter("pass: SRZjets metRel > 80    ", n_pass_SRZjets_metrel80, w);
    printCounter("pass: SRZjets 0.3<dRll<1.5   ", n_pass_SRZjets_dRll, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: CRTopmT2 OF            ", n_pass_CRTopmt2_OF, w);
    printCounter("pass: CRTopmT2 >= 1 bjet     ", n_pass_CRTopmt2_1bjet, w);
    printCounter("pass: CRTopmT2 L20/F30 veto  ", n_pass_CRTopmt2_lfveto, w);
    printCounter("pass: CRTopmT2 PtCut         ", n_pass_CRTopmt2_PtCut, w);
    printCounter("pass: CRTopmT2 mT2 > 70      ", n_pass_CRTopmt2_mt2_70, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: CRTopMet OF            ", n_pass_CRTopMet_OF, w);
    printCounter("pass: CRTopMet nB20>=1       ", n_pass_CRTopMet_1bjet, w);
    printCounter("pass: CRTopMet L20/F30 veto  ", n_pass_CRTopMet_lfveto, w);
    printCounter("pass: CRTopMet PtCut         ", n_pass_CRTopMet_PtCut, w);
    printCounter("pass: CRTopMet m(ll) < 120   ", n_pass_CRTopMet_mll120, w);
    printCounter("pass: CRTopMet pt(ll) > 80   ", n_pass_CRTopMet_ptll80, w);
    printCounter("pass: CRTopMet metRel > 80   ", n_pass_CRTopMet_metrel80, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: CRTopZjets SF          ", n_pass_CRTopZjets_SF, w);
    printCounter("pass: CRTopZjets Z Veto      ", n_pass_CRTopZjets_zv, w);
    printCounter("pass: CRTopZjets B20+L20>=2  ", n_pass_CRTopZjets_2jets, w);
    printCounter("pass: CRTopZjets B20 > 0     ", n_pass_CRTopZjets_bjet, w);
    printCounter("pass: CRTopZjets F30 = 0     ", n_pass_CRTopZjets_fveto, w);
    printCounter("pass: CRTopZjets PtCut       ", n_pass_CRTopZjets_PtCut, w);
    printCounter("pass: CRTopZjets pt(ll) > 80 ", n_pass_CRTopZjets_ptll80, w);
    printCounter("pass: CRTopZjets dRll window ", n_pass_CRTopZjets_dRll, w);
    printCounter("pass: CRTopZjets metRel > 80 ", n_pass_CRTopZjets_metrel80, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: CRWWMet OF             ", n_pass_CRWWMet_OF, w);
    printCounter("pass: CRWWMet Jet Veto       ", n_pass_CRWWMet_jv, w);
    printCounter("pass: CRWWMet PtCut          ", n_pass_CRWWMet_PtCut, w);
    printCounter("pass: CRWWMet m(ll) < 120    ", n_pass_CRWWMet_mll120, w);
    printCounter("pass: CRWWMet pt(ll) > 40    ", n_pass_CRWWMet_ptll40, w);
    printCounter("pass: CRWWMet 60<metRel<80   ", n_pass_CRWWMet_metrel_60_80, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: WWmT2 OF               ", n_pass_CRWWmt2_OF, w);
    printCounter("pass: WWmT2 Jet Veto         ", n_pass_CRWWmt2_jv, w);
    printCounter("pass: WWmT2 PtCut            ", n_pass_CRWWmt2_PtCut, w);
    printCounter("pass: WWmT2 50<mt2<90        ", n_pass_CRWWmt2_mt2_50_90, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: ZVMet Z window         ", n_pass_CRZVMet_zw, w);
    printCounter("pass: ZVMet Jet Veto         ", n_pass_CRZVMet_jv, w);
    printCounter("pass: ZVMet PtCut            ", n_pass_CRZVMet_PtCut, w);
    printCounter("pass: ZVMet pt(ll)>80        ", n_pass_CRZVMet_ptll80, w);
    printCounter("pass: ZVMet metRel > 80      ", n_pass_CRZVMet_metrel80, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: ZVmT2 Z window         ", n_pass_CRZVmt2a_zw, w);
    printCounter("pass: ZVmT2 Jet Veto         ", n_pass_CRZVmt2a_jv, w);
    printCounter("pass: ZVmT2 PtCut            ", n_pass_CRZVmt2a_PtCut, w);
    printCounter("pass: ZVmT2 mt2 > 90         ", n_pass_CRZVmt2a_mt2_90, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: ZVmT2 mt2 > 120        ", n_pass_CRZVmt2b_mt2_120, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: ZVmT2 mt2 > 150        ", n_pass_CRZVmt2c_mt2_150, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: ZVmT2 mt2 > 100        ", n_pass_CRZVmt2d_mt2_100, w);
    cout << "-----------------------------------------------------"   << endl;
    printCounter("pass: CRZXZjets SF           ", n_pass_CRZXZjets_SF, w);
    printCounter("pass: CRZXZjets Z window     ", n_pass_CRZXZjets_zw, w);
    printCounter("pass: CRZXZjets NL20>=2      ", n_pass_CRZXZjets_2ljets, w);
    printCounter("pass: CRZXZjets B20+F30=0    ", n_pass_CRZXZjets_bfveto, w);
    printCounter("pass: CRZXZjets JetPtCut     ", n_pass_CRZXZjets_JetPtCut, w);
    printCounter("pass: CRZXZjets PtCut        ", n_pass_CRZXZjets_PtCut, w);
    printCounter("pass: CRZXZjets dRll window  ", n_pass_CRZXZjets_dRll, w);
    printCounter("pass: CRZXZjets metRel > 80  ", n_pass_CRZXZjets_metrel, w);
    printCounter("pass: CRZXZjets Lower Bound  ", n_pass_CRZXZjets_lowerbound, w);
    cout << "-----------------------------------------------------"   << endl;
    cout << "Cut                   \tee\t\tmm\t\tem" << endl;
    printCounter("pass: WHSS 2lss    ", n_pass_CRWHSS2lss  , w);
    printCounter("pass: WHSS tauveto ", n_pass_CRWHSStauv  , w);
    printCounter("pass: WHSS muiso   ", n_pass_CRWHSSmuiso , w);
    printCounter("pass: WHSS eld0    ", n_pass_CRWHSSeled0 , w);
    printCounter("pass: WHSS fjveto  ", n_pass_CRWHSSnfj   , w);
    printCounter("pass: WHSS bjveto  ", n_pass_CRWHSSnbj   , w);
    printCounter("pass: WHSS nj      ", n_pass_CRWHSSnj    , w);
    printCounter("pass: WHSS 2lpt    ", n_pass_CRWHSS2lpt  , w);
    printCounter("pass: WHSS zveto   ", n_pass_CRWHSSzveto , w);
    printCounter("pass: WHSS mwwt    ", n_pass_CRWHSSmwwt  , w);
    printCounter("pass: WHSS htmin   ", n_pass_CRWHSShtmin , w);
    printCounter("pass: WHSS metrel  ", n_pass_CRWHSSmetrel, w);
    printCounter("pass: WHSS         ", n_pass_CRWHSS      , w);
  }// end loop over weight type

}

/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::printCounter(string cut, float counter[ET_N][WT_N], int weight)
{
  cout << cut;
  for(int i=0; i<ET_N-2; ++i)
    cout << "\t" << Form("%10.3f",counter[i][weight]);
  cout << endl;
}
/*--------------------------------------------------------------------------------*/
int SusySelectionMatt::getChan(const LeptonVector& leps)
{
  uint ie = 0;
  uint im = 0;
  for(uint i=0; i<leps.size(); ++i){
    if( leps.at(i)->isEle() ) ie++;
    else if( leps.at(i)->isMu() ) im++;
  }
  if( ie == 2 && im == 0 ) return Ch_ee;
  if( ie == 1 && im == 1 ) return Ch_em;
  if( ie == 0 && im == 2 ) return Ch_mm;
  cout<<"Not ee/mm/em... Number Electrons: "<<ie<<" Number Muons: "<<im<<endl;
  return Ch_N; // not in range
}
