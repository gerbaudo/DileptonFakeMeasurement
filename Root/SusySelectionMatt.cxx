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
  m_fileName("default"),
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

  //out.open(("vr1Dump_n0140/"+m_sample+".txt").c_str());
  //out.open(("vr1Dump_v2/"+m_sample+".txt").c_str());
  //out.open(("vr3Dump/"+m_sample+".txt").c_str());
  //out.open(("vr1Dump_ptWindow/"+m_sample+".txt").c_str());
  //out.open(("ZJetsDump/"+m_sample+".txt").c_str());
}

/*--------------------------------------------------------------------------------*/
// Main process loop function
/*--------------------------------------------------------------------------------*/
Bool_t SusySelectionMatt::Process(Long64_t entry)
{
  // Communicate tree entry number to SusyNtObject
  GetEntry(entry);
  clearObjects();
  m_ET = ET_Unknown;
  increment(n_readin);


  // Chain entry not the same as tree entry
  //static Long64_t chainEntry = -1;
  m_chainEntry++;
  //if( !debugEvent() ) return kTRUE;
  if(m_dbg || m_chainEntry%50000==0)
  {
    cout << "**** Processing entry " << setw(6) << m_chainEntry
         << " run " << setw(6) << nt.evt()->run
         << " event " << setw(7) << nt.evt()->event << " ****" << endl;
  }


  //cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"<<endl;
  //cout << "Run " << setw(6) << nt.evt()->run
  //<< " event " << setw(7) << nt.evt()->event << " ****" << endl;
  //((TChain*) m_tree)->GetFile()->Print();
  //return kTRUE;
  // select signal objects

  selectObjects(NtSys_NOM, false, TauID_medium);

  // Check Event
  if(!selectAnaEvent(m_signalLeptons, m_baseLeptons,true)) return kTRUE;

  //--- NEW ---//
  // Dump intersting events
  //dumpInterestingEvents(m_signalLeptons,m_signalJets2Lep, m_met);
  //return kTRUE;


  // Count SS and OS
  if(nt.evt()->isMC && !isTrueDilepton(m_signalLeptons)) return kTRUE;
  //if(nt.evt()->isMC && !(isHFLepton(m_signalLeptons[0]) || isHFLepton(m_signalLeptons[1])) )
  //return kTRUE;
  increment(n_pass_truth[m_ET], true);

  if(sameSign(m_signalLeptons))     increment(n_pass_ss[m_ET], true);

  if(oppositeSign(m_signalLeptons)) increment(n_pass_os[m_ET], true);


  // Check SR
  passSRmT2a(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passSRmT2b(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passSRmT2c(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passSRWWa(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passSRWWb(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passSRWWc(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passSRZjets(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRTopmT2(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRTopMet(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRTopZjets(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRWWMet(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRWWmT2(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRZVMet(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRZVmT2a(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRZVmT2b(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRZVmT2c(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRZVmT2d(m_signalLeptons, m_signalJets2Lep, m_met, true);
  passCRZXZjets(m_signalLeptons, m_signalJets2Lep, m_met, true);

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::Terminate()
{
  SusyNtAna::Terminate();
  if(m_dbg) cout << "SusySelectionMatt::Terminate" << endl;
  if(m_dumpCounts)
    dumpEventCounters();
}

/*--------------------------------------------------------------------------------*/
// Full event selection
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

  //--- NEW ---//
  //if( !nt.evt()->passMllForAlpgen ) return false;

  // If we are not counting (ie doing cutflow)
  // then remove all the baseline muons with
  // eta < 2.5
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
  //out<<nt.evt()->run<<" "<<nt.evt()->event<<endl;

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

  // Check signal muon if we are counting
  // Otherwise this is handled in the trigger package
  if( count && ( m_ET == ET_em || m_ET == ET_mm) )
    for(uint im=0; im<leptons.size(); ++im)
      if( leptons.at(im)->isMu() && fabs(leptons.at(im)->Eta()) > 2.4 ) return false;

  // Signal Lepotn Cut
  if( !passNLepCut(leptons) )              return false;

  // Trigger Requirement
  if( !passTrigger(baseLeps) ) return false;

  // Reject if signal taus
  if( m_signalTaus.size() != 0 )            return false;
  if(count) increment(n_pass_signalTau[m_ET]);


  // For debugging
  if(false && m_ET == ET_mm){
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
// SR mT2 a
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
// SR mT2 b
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
// SR mT2 c
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
// SR WW a
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
// SR WW b
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
// SR WW c
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
// SR Zjets
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

  // Mass requirement
  //float mjj = (*jets[0]+*jets[1]).M();
  //if( !(50 < mjj && mjj < 100) )      return false;

  // Reverse dR requirement
  //float dRll = leptons[0]->DeltaR(*leptons[1]);
  //if( (0.3 < dRll && dRll < 1.5) )              return false;

  return true;

}

/*--------------------------------------------------------------------------------*/
// Validation region SS
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
// CR WW MET
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
// CR WW mT2
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
// CR Top MET
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
// CR Top mT2
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
// CR Top Zjets
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

  //int N_L20 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  //if( N_L20 + N_B20 < 2 )                   return false;
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
// CR Top Zjets
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
// CR ZV MET
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
// CR ZV mT2
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
// Pre-mT2 region
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
bool SusySelectionMatt::passWhSS(const LeptonVector& leptons, const JetVector& jets, const Met* met)
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
  bool lsf(false), bsf(false); // compute trigw and btagw only when accepting the event
  if(sameSign(leptons)) increment(n_pass_CRWHSS2lss  [m_ET], lsf, bsf); else return false;
  DiLepEvtType ll = m_ET = getDiLepEvtType(leptons);
  bool isee(m_ET==ET_ee), isem(m_ET==ET_em||m_ET==ET_me), ismm(m_ET==ET_mm);
  float ptL0Min  = 30;
  float ptL1Min  = (ismm ? 0.0 : 20.0);
  float muIsoMax = 0.1;
  float htMin    = 200;
  float d0SMax   = (isee || isem ?   3 : FLT_MAX);
  bool applyMllZveto(isee);
  float mZ0(91.2);
  float loMllZ(applyMllZveto ? mZ0-10. : FLT_MAX);
  float hiMllZ(applyMllZveto ? mZ0+10. : FLT_MIN);
  float mtwwMin = (isee ? 150 : (isem ? 140 : (ismm ? 100 : FLT_MIN))); // todo : for now keep it simple, just one cut
  float metRelMin = (isee ? 50 : (isem ? 50 : (ismm ? FLT_MIN : FLT_MIN))); // for now simple

  if(m_signalTaus.size()==0)               increment(n_pass_CRWHSStauv  [m_ET], lsf, bsf); else return false;
  if(numberOfFJets(jets)==0)                          increment(n_pass_CRWHSSnfj   [m_ET], lsf, bsf); else  return false;
  if(numberOfCBJets(jets)==0)                         increment(n_pass_CRWHSSnbj   [m_ET], lsf, bsf); else  return false;
  if(numberOfCLJets(jets)>0)                          increment(n_pass_CRWHSSnj    [m_ET], lsf, bsf); else  return false;
  if(susy::pass2LepPt    (leptons, ptL0Min, ptL1Min)) increment(n_pass_CRWHSS2lpt  [m_ET], lsf, bsf); else  return false;
  if(susy::passZllVeto   (leptons, loMllZ, hiMllZ))   increment(n_pass_CRWHSSzveto [m_ET], lsf, bsf); else  return false;
  if(susy::passMtLlMetMin(leptons, met, mtwwMin))     increment(n_pass_CRWHSSmwwt  [m_ET], lsf, bsf); else  return false;
  if(susy::passHtMin     (leptons, jets, met, htMin)) increment(n_pass_CRWHSShtmin [m_ET], lsf, bsf); else  return false;
  if(getMetRel(met,leptons,jets)>metRelMin)           increment(n_pass_CRWHSSmetrel[m_ET], lsf, bsf); else  return false;
  lsf = bsf = true;
  increment(n_pass_CRWHSS[m_ET], lsf, bsf);
  return true;
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

  //cout<<"NLeptons: "<<leptons.size()<<endl;
  //for(uint i=0; i<leptons.size(); ++i)
  //cout<<"\tElectron: "<<leptons.at(i)->isEle()<<" "<<leptons.at(i)->trigFlags<<endl;

  if(leptons.size() != 2)
    return false;

  // Trigger Check -- muon only events
  //if( !(leptons[0]->isMu() && leptons[1]->isMu()) ) return false;
  //if( leptons[0]->Pt() < 14 ) return false;
  //if( leptons[1]->Pt() < 14 ) return false;

  //bool passEvtTrig = nt.evt()->passTrig( TRIG_2mu13 );
  //bool passTrigMatch = leptons[0]->matchTrig( TRIG_2mu13 ) && leptons[1]->matchTrig( TRIG_2mu13 );

  bool passEvtTrig   = m_trigObj->passDilEvtTrig(leptons, m_met->Et, nt.evt());
  bool passTrigMatch = m_trigObj->passDilTrigMatch(leptons, m_met->Et, nt.evt());

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
bool SusySelectionMatt::sameFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return (leptons.at(0)->isMu() == leptons.at(1)->isMu());
  //return (leptons.at(0)->isEle() && leptons.at(1)->isEle());
  //return (leptons.at(0)->isMu() && leptons.at(1)->isMu());
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
// Signal region cuts
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passJetVeto(const JetVector& jets)
{

  // Require no light, b, or forward jets
  int N_L25 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);

  return (N_L25 + N_B20 + N_F30 == 0);

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
bool SusySelectionMatt::passge2Jet(const JetVector& jets)
{

  // N_L25 >=2 N_B20 = N_F30 = 0
  int N_L25 = numberOfCLJets(jets);
  int N_B20 = numberOfCBJets(jets);
  int N_F30 = numberOfFJets(jets);

  return (N_L25 >=2 && N_B20 + N_F30 == 0);

  /*
  // Count jets above 30 GeV
  int njet = 0;
  for(uint j=0; j<jets.size(); ++j)
    if( jets.at(j)->Pt() > 30 )
      njet++;

  return (njet >= 2);
  */
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

  // necessary variables
  TLorentzVector metlv = met->lv();
  TLorentzVector l0    = *leptons.at(0);
  TLorentzVector l1    = *leptons.at(1);

  double pTMiss[3] = {0.0, metlv.Px(), metlv.Py()};
  double pA[3]     = {0.0, l0.Px(), l0.Py()};
  double pB[3]     = {0.0, l1.Px(), l1.Py()};

  // Create Mt2 object
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pA,pB,pTMiss);
  mt2_event.set_mn(0); // LSP mass = 0 is Generic

  return mt2_event.get_mt2();
}

//-----------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/InDetTrackingPerformanceGuidelines
// MC beam spot size sigma_z = 66mm
// Data beam spot size sigma_z = 47mm
//Retunr the number of matching Npv in data given a MC
// I copied this from Anyes
float SusySelectionMatt::GetNVertexBsCorrected(float nRecoVtx)
{

  TGraph *g_nvtx_nreco_bs66mm = new TGraph(41);
  g_nvtx_nreco_bs66mm->SetName("g_shadowedAverage_bs66mm");
  g_nvtx_nreco_bs66mm->SetPoint(0, 0, 0);
  g_nvtx_nreco_bs66mm->SetPoint(1, 1, 1);
  g_nvtx_nreco_bs66mm->SetPoint(2, 2, 1.97943);
  g_nvtx_nreco_bs66mm->SetPoint(3, 3, 2.93912);
  g_nvtx_nreco_bs66mm->SetPoint(4, 4, 3.87986);
  g_nvtx_nreco_bs66mm->SetPoint(5, 5, 4.8024);
  g_nvtx_nreco_bs66mm->SetPoint(6, 6, 5.70743);
  g_nvtx_nreco_bs66mm->SetPoint(7, 7, 6.59561);
  g_nvtx_nreco_bs66mm->SetPoint(8, 8, 7.46756);
  g_nvtx_nreco_bs66mm->SetPoint(9, 9, 8.32386);
  g_nvtx_nreco_bs66mm->SetPoint(10, 10, 9.16509);
  g_nvtx_nreco_bs66mm->SetPoint(11, 11, 9.99175);
  g_nvtx_nreco_bs66mm->SetPoint(12, 12, 10.8044);
  g_nvtx_nreco_bs66mm->SetPoint(13, 13, 11.6034);
  g_nvtx_nreco_bs66mm->SetPoint(14, 14, 12.3893);
  g_nvtx_nreco_bs66mm->SetPoint(15, 15, 13.1624);
  g_nvtx_nreco_bs66mm->SetPoint(16, 16, 13.9233);
  g_nvtx_nreco_bs66mm->SetPoint(17, 17, 14.6723);
  g_nvtx_nreco_bs66mm->SetPoint(18, 18, 15.4097);
  g_nvtx_nreco_bs66mm->SetPoint(19, 19, 16.136);
  g_nvtx_nreco_bs66mm->SetPoint(20, 20, 16.8514);
  g_nvtx_nreco_bs66mm->SetPoint(21, 21, 17.5562);
  g_nvtx_nreco_bs66mm->SetPoint(22, 22, 18.2509);
  g_nvtx_nreco_bs66mm->SetPoint(23, 23, 18.9356);
  g_nvtx_nreco_bs66mm->SetPoint(24, 24, 19.6106);
  g_nvtx_nreco_bs66mm->SetPoint(25, 25, 20.2763);
  g_nvtx_nreco_bs66mm->SetPoint(26, 26, 20.9329);
  g_nvtx_nreco_bs66mm->SetPoint(27, 27, 21.5806);
  g_nvtx_nreco_bs66mm->SetPoint(28, 28, 22.2197);
  g_nvtx_nreco_bs66mm->SetPoint(29, 29, 22.8504);
  g_nvtx_nreco_bs66mm->SetPoint(30, 30, 23.4729);
  g_nvtx_nreco_bs66mm->SetPoint(31, 31, 24.0874);
  g_nvtx_nreco_bs66mm->SetPoint(32, 32, 24.6941);
  g_nvtx_nreco_bs66mm->SetPoint(33, 33, 25.2933);
  g_nvtx_nreco_bs66mm->SetPoint(34, 34, 25.8851);
  g_nvtx_nreco_bs66mm->SetPoint(35, 35, 26.4697);
  g_nvtx_nreco_bs66mm->SetPoint(36, 36, 27.0473);
  g_nvtx_nreco_bs66mm->SetPoint(37, 37, 27.6179);
  g_nvtx_nreco_bs66mm->SetPoint(38, 38, 28.1819);
  g_nvtx_nreco_bs66mm->SetPoint(39, 39, 28.7394);
  g_nvtx_nreco_bs66mm->SetPoint(40, 40, 29.2904);

  TGraph *g_nvtx_nreco_bs47mm = new TGraph(41);
  g_nvtx_nreco_bs47mm->SetName("g_shadowedAverage_bs47mm");
  g_nvtx_nreco_bs47mm->SetPoint(0, 0, 0);
  g_nvtx_nreco_bs47mm->SetPoint(1, 1, 1);
  g_nvtx_nreco_bs47mm->SetPoint(2, 2, 1.97111);
  g_nvtx_nreco_bs47mm->SetPoint(3, 3, 2.91499);
  g_nvtx_nreco_bs47mm->SetPoint(4, 4, 3.83313);
  g_nvtx_nreco_bs47mm->SetPoint(5, 5, 4.72692);
  g_nvtx_nreco_bs47mm->SetPoint(6, 6, 5.59763);
  g_nvtx_nreco_bs47mm->SetPoint(7, 7, 6.44645);
  g_nvtx_nreco_bs47mm->SetPoint(8, 8, 7.27445);
  g_nvtx_nreco_bs47mm->SetPoint(9, 9, 8.08265);
  g_nvtx_nreco_bs47mm->SetPoint(10, 10, 8.87199);
  g_nvtx_nreco_bs47mm->SetPoint(11, 11, 9.64332);
  g_nvtx_nreco_bs47mm->SetPoint(12, 12, 10.3975);
  g_nvtx_nreco_bs47mm->SetPoint(13, 13, 11.1352);
  g_nvtx_nreco_bs47mm->SetPoint(14, 14, 11.8572);
  g_nvtx_nreco_bs47mm->SetPoint(15, 15, 12.5641);
  g_nvtx_nreco_bs47mm->SetPoint(16, 16, 13.2567);
  g_nvtx_nreco_bs47mm->SetPoint(17, 17, 13.9353);
  g_nvtx_nreco_bs47mm->SetPoint(18, 18, 14.6007);
  g_nvtx_nreco_bs47mm->SetPoint(19, 19, 15.2532);
  g_nvtx_nreco_bs47mm->SetPoint(20, 20, 15.8935);
  g_nvtx_nreco_bs47mm->SetPoint(21, 21, 16.5219);
  g_nvtx_nreco_bs47mm->SetPoint(22, 22, 17.1389);
  g_nvtx_nreco_bs47mm->SetPoint(23, 23, 17.7449);
  g_nvtx_nreco_bs47mm->SetPoint(24, 24, 18.3404);
  g_nvtx_nreco_bs47mm->SetPoint(25, 25, 18.9255);
  g_nvtx_nreco_bs47mm->SetPoint(26, 26, 19.5008);
  g_nvtx_nreco_bs47mm->SetPoint(27, 27, 20.0665);
  g_nvtx_nreco_bs47mm->SetPoint(28, 28, 20.623);
  g_nvtx_nreco_bs47mm->SetPoint(29, 29, 21.1705);
  g_nvtx_nreco_bs47mm->SetPoint(30, 30, 21.7094);
  g_nvtx_nreco_bs47mm->SetPoint(31, 31, 22.2399);
  g_nvtx_nreco_bs47mm->SetPoint(32, 32, 22.7622);
  g_nvtx_nreco_bs47mm->SetPoint(33, 33, 23.2767);
  g_nvtx_nreco_bs47mm->SetPoint(34, 34, 23.7835);
  g_nvtx_nreco_bs47mm->SetPoint(35, 35, 24.2829);
  g_nvtx_nreco_bs47mm->SetPoint(36, 36, 24.7751);
  g_nvtx_nreco_bs47mm->SetPoint(37, 37, 25.2604);
  g_nvtx_nreco_bs47mm->SetPoint(38, 38, 25.7388);
  g_nvtx_nreco_bs47mm->SetPoint(39, 39, 26.2106);
  g_nvtx_nreco_bs47mm->SetPoint(40, 40, 26.6759);


  //get corresponding NReconstructible (points are already sorted in X, monotonic in Y)
  float nRecon=-1;
  if (nRecoVtx < g_nvtx_nreco_bs66mm->GetY()[0]) {
    //std::cout << "NVtx_bs_correction.C: WARNING: Requested nVertex outside the expected range: " << nRecoVtx << std::endl;
    return nRecoVtx; //do not correct
  }
  if (nRecoVtx > g_nvtx_nreco_bs66mm->GetY()[g_nvtx_nreco_bs66mm->GetN()-1]) {
    //std::cout << "NVtx_bs_correction.C: WARNING: Requested nVertex outside the expected range: " << nRecoVtx << std::endl;
    return nRecoVtx; //do not correct
  }
  for (int i=1; i < g_nvtx_nreco_bs66mm->GetN(); i++) {
    if (nRecoVtx < g_nvtx_nreco_bs66mm->GetY()[i]) {
      //linear interpolation
      nRecon = g_nvtx_nreco_bs66mm->GetX()[i-1]+(nRecoVtx - (g_nvtx_nreco_bs66mm->GetY()[i-1])) *
	(g_nvtx_nreco_bs66mm->GetX()[i] - g_nvtx_nreco_bs66mm->GetX()[i-1]) /
	(g_nvtx_nreco_bs66mm->GetY()[i] - g_nvtx_nreco_bs66mm->GetY()[i-1]);

      break;
    }
  }

  //now return corresponding reconstructed vertices for bs=47mm

  float val = g_nvtx_nreco_bs47mm->Eval(nRecon);

  delete g_nvtx_nreco_bs47mm;
  delete g_nvtx_nreco_bs66mm;

  return val;

}


/*--------------------------------------------------------------------------------*/
// Check MC Lepton
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

  //uint origin = lep->mcOrigin;
  //return origin == 25 ||origin == 26 || origin == 27 || origin == 28 ||
  //origin == 29 || origin == 32 || origin == 33;



}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isLFLepton(const Lepton* lep)
{

  return (lep->truthType == RecoTruthMatch::LF);

  // Steve's way:
  //bool isChargeFlip = lep->isEle() ? ((Electron*) lep)->isChargeFlip : false;
  //return isFakeLepton(lep) && !isConvLepton(lep) && !isHFLepton(lep) && !isChargeFlip;

  // 2-lep way:
  //uint origin = lep->mcOrigin;
  //return origin == 0 || origin == 23 || origin == 24 || origin == 30 ||
  //origin == 31 || origin == 34 || origin == 35;


}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isQCDLepton(const Lepton* lep)
{

  return isHFLepton(lep) || isLFLepton(lep);

}
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isTrueDilepton(const LeptonVector &leptons)
{

  //if( leptons.size() !=2 ) return false;
  //  bool l0real = leptons[0]->truthType == RecoTruthMatch::PROMPT;
  //bool l1real = leptons[1]->truthType == RecoTruthMatch::PROMPT;
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
bool SusySelectionMatt::isFakeDilepton(const LeptonVector &leptons)
{

  // Maybe not 100% kosher, but I just want to make sure I am not
  // double counting, so I require dilepton events to be real
  if( leptons.size() != 2 ) return false;

  //cout<<"Lepton 0: "<<isFakeLepton(leptons[0])<<" truth: "<<leptons[0]->truthType<<" type: "<<leptons[0]->mcType<<" origin: "<<leptons[1]->mcOrigin<<endl;
  //cout<<"Lepton 1: "<<isFakeLepton(leptons[1])<<" truth: "<<leptons[1]->truthType<<" type: "<<leptons[1]->mcType<<" origin: "<<leptons[1]->mcOrigin<<endl;


  bool l0_fake = isFakeLepton(leptons[0]);
  bool l1_fake = isFakeLepton(leptons[1]);
  bool l0_cf = leptons[0]->isEle() ? ((Electron*) leptons[0])->isChargeFlip : false;
  bool l1_cf = leptons[1]->isEle() ? ((Electron*) leptons[1])->isChargeFlip : false;

  return (l0_fake || l1_fake) && !l0_cf && !l1_cf; // ignoring charge flip
  //return (l0_real || l0_cf) && (l1_real || l1_cf);

}

/*--------------------------------------------------------------------------------*/
// Get Event weight
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
  //else weight = getEventWeightFixed(nt.evt()->mcChannel, LUMI_A_L);
  else weight = getEventWeight(LUMI_A_L, true);
  //else weight = getEventWeight(LUMI_A_B14, true);

  //if(m_do1fb) weight = getEventWeightFixed(nt.evt()->mcChannel, LUMI_A_B3);
  //else if(m_doAD)  weight = getEventWeightFixed(nt.evt()->mcChannel,LUMI_A_D);
  //else weight = getEventWeightFixed(nt.evt()->mcChannel,LUMI_A_E);

  // bbbar/ccbar scale factor
  /*
  uint chNum = nt.evt()->mcChannel;
  if(chNum == 129136 || chNum == 147668){
    if(leptons[0]->isMu() && leptons[1]->isMu()) weight *= 0.706074;
    else weight = 0;
  }
  if(chNum == 129135 || chNum == 147667){
    if(leptons[0]->isEle() && leptons[1]->isEle()) weight *= 0.706074;
    else weight = 0;
  }
  */

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
float SusySelectionMatt::getBTagWeight(const Event* evt)
{
  JetVector tempJets;
  for(uint ij=0; ij<m_baseJets.size(); ++ij){
    Jet* jet = m_baseJets.at(ij);
    //if( !(jet->Pt() > 20 && fabs(jet->Eta()) < 2.4) ) continue;
    if( !(jet->Pt() > 20 && fabs(jet->detEta) < 2.4) ) continue;
    tempJets.push_back(jet);
  }

  //cout<<"Nominal: "
  //<<bTagSF(evt, tempJets, nt.evt()->mcChannel, BTag_NOM)
  //<<" Up: "
  //<<bTagSF(evt, tempJets, nt.evt()->mcChannel, BTag_BJet_UP)
  //<<" Down: "
  //<<bTagSF(evt, tempJets, nt.evt()->mcChannel, BTag_BJet_DN)
  //<<endl;

  //return bTagSF(evt, tempJets, true, "MV1", "0_3511", MV1_80, BTag_NOM);
  return bTagSF(evt, tempJets, evt->mcChannel, BTag_NOM);
}

/*--------------------------------------------------------------------------------*/
// Photon methods
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::passCheckMC(int mcRunNumber, float pt)
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
float SusySelectionMatt::getPhotonXS(int mcRunNumber)
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
void SusySelectionMatt::dumpEventCounters()
{


  //string v_ET[] = {"ee","mm","em","me"};
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
// Print for given di-lepton category
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::printCounter(string cut, float counter[ET_N][WT_N], int weight)
{
  cout << cut;
  for(int i=0; i<ET_N-2; ++i)
    cout << "\t" << Form("%10.3f",counter[i][weight]);
  cout << endl;
}

/*--------------------------------------------------------------------------------*/
// Dump Pre objects
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::dumpPreObjects()
{
  ElectronVector preElectrons = getPreElectrons(&nt, NtSys_NOM);
  MuonVector preMuons = getPreMuons(&nt, NtSys_NOM);
  JetVector preJets = getPreJets(&nt, NtSys_NOM);

  out << "Pre Electrons: "<<preElectrons.size()<<endl;
  for(uint ie=0; ie<preElectrons.size(); ++ie){
    preElectrons[ie]->print();
    //out<<"Eta: "<<preElectrons[ie]->Eta()<<" Phi: "<<preElectrons[ie]->Phi()<<endl;
  }

  out<<endl;
  out << "Pre Muons: "<<preMuons.size() << endl;
  for(uint im=0; im<preMuons.size(); ++im){
    preMuons[im]->print();
    out<<"ptcone: "<<preMuons[im]->ptcone30<<" corrected: "<<muPtConeCorr(preMuons[im],m_baseElectrons, m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC)<<endl;
    out<<"ptcone/pt: "<<muPtConeCorr(preMuons[im],m_baseElectrons, m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC)/preMuons[im]->Pt()<<endl;
    out<<"d0: "<<preMuons[im]->d0Unbiased<<" Err d0 unb: "<<preMuons[im]->errD0Unbiased<<endl;
    out<<"z0: "<<preMuons[im]->z0Unbiased<<" Err z0 unb: "<<preMuons[im]->errZ0Unbiased<<endl;
    out<<"d0sig: "<<preMuons[im]->d0Sig()<<" z0*sin(theta): "<<preMuons[im]->z0SinTheta()<<endl;

    //out<<"Eta: "<<preMuons[im]->Eta()<<" Phi: "<<preMuons[im]->Phi()<<endl;
  }

  out << endl;
  out << "Pre Jets: "<<preJets.size() << endl;
  for(uint ij=0; ij<preJets.size(); ++ij){
    preJets[ij]->print();
    //out<<"Eta: "<<preJets[ij]->Eta()<<" Phi: "<<preJets[ij]->Phi()<<endl;
  }

}
/*--------------------------------------------------------------------------------*/
// Dump Jets with more information
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::dumpJets()
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
void SusySelectionMatt::checkSys()
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
bool SusySelectionMatt::debugEvent()
{
  uint run = nt.evt()->run;
  uint evt = nt.evt()->event;

  /*
  if(run==204265 && evt ==99597041) return true;
  if(run==203195 && evt ==41438450) return true;
  if(run==204240 && evt ==33857176) return true;
  if(run==203456 && evt ==14297388) return true;
  if(run==208970 && evt ==77422707) return true;
  if(run==208354 && evt ==28241368) return true;
  if(run==207620 && evt ==93663162) return true;
  if(run==207620 && evt ==94247317) return true;
  if(run==208662 && evt ==136103149) return true;
  if(run==208258 && evt ==2678576) return true;
  if(run==209183 && evt ==187248365) return true;
  if(run==209109 && evt ==78161028) return true;
  if(run==209629 && evt ==185790707) return true;
  if(run==209608 && evt ==45653361) return true;
  */
  //if(evt == 1623445 ) return true;
  //if(evt == 84903) return true;
  //if(evt == 801176) return true;
  //if(evt == 840519) return true;
  //if(evt == 378326) return true;
  //if(evt == 848393) return true;
  //if(evt == 97552) return true;
  //if(evt == 378326) return true;

  if(run==195847 && evt==1094861)return true;
  if(run==195847 && evt==1108285)return true;
  if(run==195847 && evt==1322916)return true;
  if(run==195847 && evt==1833312)return true;
  if(run==195847 && evt==2578769)return true;
  if(run==195847 && evt==670066)return true;
  if(run==195847 && evt==739943)return true;
  if(run==195847 && evt==782050)return true;


  return false;
}

/*--------------------------------------------------------------------------------*/
// Dump interesting events
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::dumpInterestingEvents(const LeptonVector& leptons,
					  const JetVector& jets,
					  const Met* met)
{

  //if( !passBR2(m_signalLeptons, m_signalJets2Lep, m_met) ) return;
  // Look for 1 jet, SS muons, 90-120 GeV
  //if(leptons.size() !=2) return;
  if( !passSRZjets(leptons, jets, met) ) return;

  float mll = Mll(leptons[0],leptons[1]);
  //if( !passVRSS(leptons, jets, met) ) return;
  //if( !passVR3(leptons, jets, met) ) return;
  if( !(leptons[0]->isMu() && leptons[1]->isMu()) ) return;
  //if( leptons[0]->isEle() && leptons[1]->isEle()) return;

  float metRel = getMetRel(met,leptons,jets);
  if( metRel < 100 ) return;

  //float pt0 = leptons[0]->Pt();
  //if( !(20 < pt0 && pt0 < 30) ) return;
  //if( !(50 < pt0 && pt0 < 60) ) return;
  /*
  if(!sameSign(leptons)) return;
  if(!(90 < mll && mll < 120)) return;
  if(jets.size() != 1) return;
  if( met->Et < 40 ) return;

  if( !(leptons[0]->isMu() && leptons[1]->isMu()) ) return;
  */
  // Now dump the interesting variables
  const Lepton* l0 = leptons[0];
  const Lepton* l1 = leptons[1];
  const Jet*    j  = jets.size() > 0 ? jets[0] : NULL;



  /*out<<"Run "<<nt.evt()->run
     <<" Event "<<nt.evt()->event
     <<" mll "<<Mll(l0,l1)
     <<" l0pt "<<l0->Pt()
     <<" l0eta "<<l0->Eta()
     <<" l1pt "<<l1->Pt()
     <<" l1eta "<<l1->Eta()
     <<" l0q "<<l0->q
     <<" l1q "<<l1->q
     <<" njets "<<jets.size()
     <<" met "<<met->Et
     <<" metrel "<<metRel
     <<endl;*/


  out<<"-----------------------------------------------"<<endl;
  out<<"Run: "<<nt.evt()->run<<" Event: "<<nt.evt()->event<<endl;
  out<<"l0: Pt "<<l0->Pt()<<" Eta "<<l0->Eta()<<" Phi "<<l0->Phi()<<" isEle: "<<l0->isEle()<<" d0Sig: "<<l0->d0Sig(true)<<" biased: "<<l0->d0Sig(false)<<" z0sintheta: "<<l0->z0SinTheta(true)<<" biased: "<<l0->z0SinTheta(false)<<endl;
  out<<"l1: Pt "<<l1->Pt()<<" Eta "<<l1->Eta()<<" Phi "<<l1->Phi()<<" isEle: "<<l1->isEle()<<" d0Sig: "<<l1->d0Sig(true)<<" biased: "<<l1->d0Sig(false)<<" z0sintheta: "<<l1->z0SinTheta(true)<<" biased: "<<l1->z0SinTheta(false)<<endl;

  // Isolation
  if( l0->isMu() ){
    out<<" l0 ptcone30: "<<((Muon*)l0)->ptcone30ElStyle<<" relative: "<<(((Muon*)l0)->ptcone30ElStyle / l0->Pt())
       <<" l0 etcone30: "<<((Muon*)l0)->etcone30<<" relative: "<<(((Muon*)l0)->etcone30 / l0->Pt())<<endl;
    out<<" l0 isCombined: "<<((Muon*)l0)->isCombined<<" idTrackQ: "<<((Muon*)l0)->idTrackQ<<" msTrackQ: "<<((Muon*)l0)->msTrackQ<<endl;
  }
  if( l1->isMu() ){
    out<<" l1 ptcone30: "<<((Muon*)l1)->ptcone30ElStyle<<" relative: "<<(((Muon*)l1)->ptcone30ElStyle / l1->Pt())
       <<" l1 etcone30: "<<((Muon*)l1)->etcone30<<" relative: "<<(((Muon*)l1)->etcone30 / l1->Pt())<<endl;
    out<<" l1 isCombined: "<<((Muon*)l1)->isCombined<<" idTrackQ: "<<((Muon*)l1)->idTrackQ<<" msTrackQ: "<<((Muon*)l1)->msTrackQ<<endl;
  }

  out<<"Njets: "<<jets.size()<<endl;
  for(uint ij=0; ij<jets.size(); ++ij)
    out<<"Jet: Pt "<<jets.at(ij)->Pt()
       <<" Eta "<<jets.at(ij)->Eta()
       <<" Phi "<<jets.at(ij)->Phi()
       <<" jvf: "<<jets.at(ij)->jvf
       <<" mv1: "<<jets.at(ij)->mv1
       <<" isCentralL: "<<isCentralLightJet(jets.at(ij))
       <<" isCentralB: "<<isCentralBJet(jets.at(ij))
       <<" isForward: "<<isForwardJet(jets.at(ij))<<endl;
  out<<"Mll: "<<mll<<" Met "<<met->Et<<" metRel "<<metRel<<" phi: "<<met->phi<<endl;

  if(j) out<<"mllj "<<(*l0 + *l1 + *j).M()<<endl;

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

  LeptonVector preLeptons;
  ElectronVector preEl = getPreElectrons(&nt, NtSys_NOM);
  MuonVector preMu = getPreMuons(&nt, NtSys_NOM);
  TauVector preTau = getPreTaus(&nt, NtSys_NOM);
  buildLeptons(preLeptons, preEl, preMu);
  out<<"-------- Pre Objects ----------"<<endl;
  out<<"Pre leptons: "<<preLeptons.size()<<endl;
  if(preLeptons.size() > 2){
    for(uint iL=0; iL < preLeptons.size(); ++iL){
      const Lepton* lep = preLeptons.at(iL);
      if( lep == l0 || lep == l1) continue;
      out<<"Lep: "<<iL<<" Pt "<<lep->Pt()<<" Eta: "<<lep->Eta()<<" Phi: "<<lep->Phi()
	 <<" isEle: "<<lep->isEle()<<endl;
      out<<"Ml0l"<<iL<<" = "<<Mll(l0,lep)<<endl;
      out<<"Ml1l"<<iL<<" = "<<Mll(l1,lep)<<endl;

    }
  }
  out<<"-------- Pre Taus ----------"<<endl;
  out<<"Pre taus: "<<preTau.size()<<endl;
  for(uint it=0; it<preTau.size(); ++it){
    const Tau* tau = preTau.at(it);
    out<<" tau"<<it<<" Pt: "<<tau->Pt()<<" Eta: "<<tau->Eta()<<" Phi: "<<tau->Phi()
       <<" ntrack "<<tau->nTrack<<" eleBDT: "<<tau->eleBDT<<" jetBDT: "<<tau->jetBDT
       <<" muVeto: "<<tau->muonVeto<<endl;
  }


}
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::dumpTrigFlag(uint flag)
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
void SusySelectionMatt::printLep(const Lepton* lep)
{

  out<<"\tLepton is Muon: "<<lep->isMu()<<" Charge: "<<lep->q<<endl;
  if(lep->isMu())out<<"\tCombined: "<<((Muon*)lep)->isCombined<<endl;
  out<<"\tPt: "<<lep->Pt()<<" Eta: "<<lep->Eta()<<" Phi: "<<lep->Phi()<<endl;
  out<<"\tptcone20: "<<lep->ptcone20<<" ptcone30: "<<lep->ptcone30<<endl;
  out<<"\td0: "<<lep->d0<<" d0err: "<<lep->errD0<<endl;

}
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::printJet(const Jet* jet)
{

  out<<"\tPt: "<<jet->Pt()<<" Eta: "<<jet->Eta()<<" Phi: "<<jet->Phi()<<endl;
  out<<"\tjvf: "<<jet->jvf<<" sv0: "<<jet->sv0
     <<" combNN: "<<jet->combNN<<" mv1: "<<jet->mv1<<endl;

}
/*--------------------------------------------------------------------------------*/
// Get lepton channel
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

/*--------------------------------------------------------------------------------*/
// Test method to modify baseline definition
/*--------------------------------------------------------------------------------*/
bool SusySelectionMatt::isBaselineLepton(const LeptonVector &leptons)
{

  // Testing if I put d0 and z0 in the denominator
  for(uint il=0; il<leptons.size(); ++il){

    float d0cut = leptons[il]->isEle() ? ELECTRON_D0SIG_CUT : MUON_D0SIG_CUT;
    float z0cut = leptons[il]->isEle() ? ELECTRON_Z0_SINTHETA_CUT : MUON_Z0_SINTHETA_CUT;

    if( fabs(leptons[il]->d0Sig(true)) >= d0cut )      return false;
    if( fabs(leptons[il]->z0SinTheta(true)) >= z0cut ) return false;

    //if( leptons[il]->isEle() && !((Electron*) leptons[il])->tightPP) return false;

  }

  return true;


}

/*--------------------------------------------------------------------------------*/
// Testing method for extending the loose definition
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::selectFakeObjects()
{

  // Want to modify slightly the baseline lepton definition to be used
  // in the fake estimate to try to increase the statistics.  Will
  // naively steal cuts from Z+b analysis

  clearObjects();

  // Manually re-write the getBaselineObjects
  m_baseElectrons = getPreElectrons(&nt, NtSys_NOM);
  m_baseMuons     = getPreMuons(&nt, NtSys_NOM);
  m_baseJets      = getPreJets(&nt, NtSys_NOM);
  m_baseTaus      = getPreTaus(&nt, NtSys_NOM);

  // Overlap removal:

  // Remove electrons from electrons
  e_e_overlap(m_baseElectrons, E_E_DR);

  // Remove jets from electrons
  e_j_overlap(m_baseElectrons, m_baseJets, J_E_DR, true);

  // Remove taus from electrons
  t_e_overlap(m_baseTaus, m_baseElectrons, T_E_DR);

  // Remove taus from muons
  t_m_overlap(m_baseTaus, m_baseMuons, T_M_DR);

  // Remove electrons from jets
  e_j_overlap(m_baseElectrons, m_baseJets, E_J_DR, false);
  // Remove muons from jets
  // Write this one to allow more muons
  if( !(m_baseMuons.size()==0 || m_baseJets.size()==0) ){

    for(int im=m_baseMuons.size()-1; im>=0; im--){
      const Muon* mu = m_baseMuons.at(im);
      for(int ij=m_baseJets.size()-1; ij>=0; ij--){
        const Jet* j = m_baseJets.at(ij);

        float dR = mu->DeltaR(*j);

        if(dR > 0.4) continue;
        if( 0.2 < dR && mu->ptcone30ElStyle/mu->Pt() < 0.3) continue;

        m_baseMuons.erase( m_baseMuons.begin() + im );
        break;

      }// end loop over jets
    }// end loop over muons

  }// End if have muons and jets

  // Remove electrons and muons that overlap
  e_m_overlap(m_baseElectrons, m_baseMuons, E_M_DR);

  // Remove muons from muons
  m_m_overlap(m_baseMuons, M_M_DR);

  // Remove jets from taus
  t_j_overlap(m_baseTaus, m_baseJets, J_T_DR, true);

  // Get the signal objects
  getSignalObjects(m_baseElectrons, m_baseMuons, m_baseTaus, m_baseJets,
                   m_signalElectrons, m_signalMuons, m_signalTaus, m_signalJets,
                   m_signalJets2Lep, nt.evt()->nVtx, nt.evt()->isMC, false);

  m_signalTaus = m_mediumTaus; // We use medium taus

  // Additionaly, the muons cannot overlap with a jet...
  MuonVector temp = m_signalMuons;
  for(uint im=0; im<temp.size(); ++im)
    if( isMySignalMuon(temp.at(im)) ) m_signalMuons.push_back( temp.at(im) );


  // Get Met
  m_met = getMet(&nt, NtSys_NOM);

  // Build lepton vectors
  buildLeptons(m_baseLeptons, m_baseElectrons, m_baseMuons);
  buildLeptons(m_signalLeptons, m_signalElectrons, m_signalMuons);

  std::sort(m_baseLeptons.begin(), m_baseLeptons.end(), comparePt);
  std::sort(m_signalLeptons.begin(), m_signalLeptons.end(), comparePt);
  std::sort(m_baseJets.begin(), m_baseJets.end(), comparePt);
  std::sort(m_signalJets2Lep.begin(), m_signalJets2Lep.end(), comparePt);

}

bool SusySelectionMatt::isMySignalMuon(const Muon* mu)
{

  if( !isSignalMuon(mu, m_baseElectrons, m_baseMuons,
                    nt.evt()->nVtx, nt.evt()->isMC, false)
      )
    return false;

  for( uint ij=0; ij<m_baseJets.size(); ++ij){
    const Jet* j = m_baseJets.at(ij);

    if( j->DeltaR(*mu) < 0.4 ) return false;

  }

  return true;

}

bool SusySelectionMatt::isMySignalLepton(const Lepton* lep)
{

  if( lep->isEle() )
    return isSignalLepton(lep,
			  m_baseElectrons,
			  m_baseMuons,
			  nt.evt()->nVtx,
			  nt.evt()->isMC,
			  false);

  return isMySignalMuon((Muon*) lep);

}
