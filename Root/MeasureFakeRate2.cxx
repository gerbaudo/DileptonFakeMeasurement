////////////////////////////////////////////////////
// Will house the control regions used to measure //
// the fake rates used in 2-lep analysis.         //
////////////////////////////////////////////////////

#include "SusyTest0/MeasureFakeRate2.h"
#include "SusyTest0/SelectionRegions.h"
#include "SusyTest0/criteria.h"

float Htbins[] = {1,3,6,9,12};
int nHtbins = 4;


using namespace susywh; // pull in SelectionRegions
const int controlRegions[] = {
  kCR_Real, kCR_SideLow, kCR_SideHigh, kCR_HF, kCR_HF_high,
  kCR_LFZjet, kCR_LFWjet, kCR_Conv, kCR_CFlip,
  kCR_MCHeavy, kCR_MCLight, kCR_MCConv, kCR_MCQCD,
  kCR_MCALL, kCR_MCReal, kCR_MCNone
};
const size_t nControlRegions = sizeof(controlRegions)/sizeof(int);
const int signalRegions[] = {
  kCR_SRmT2a,
  kCR_SRmT2b,
  kCR_SRmT2c,
  kCR_SRWWa,
  kCR_SRWWb,
  kCR_SRWWc,
  kCR_SRZjets,
  kCR_VRSS,
  kCR_CRWWMet,
  kCR_CRWWmT2,
  kCR_CRTopMet,
  kCR_CRTopmT2,
  kCR_CRZVMet,
  kCR_CRZVmT2_90,
  kCR_CRZVmT2_120,
  kCR_CRZVmT2_150,
  kCR_CRZVmT2_100,
  kCR_CRTopZjets,
  kCR_CRZXZjets,
  kCR_SSInc,
  kCR_PremT2,
  kCR_SRWHSS,
};
const size_t nSignalRegions = sizeof(signalRegions)/sizeof(int);

/*--------------------------------------------------------------------------------*/
// Fake Rate Constructor
/*--------------------------------------------------------------------------------*/
MeasureFakeRate2::MeasureFakeRate2() :
  CR_N(nControlRegions + nSignalRegions + 1),
  m_controlRegions(controlRegions, controlRegions + nControlRegions),
  m_signalRegions (signalRegions,  signalRegions  + nSignalRegions),
  m_outFile(NULL),
  m_evtWeight(1.),
  m_metRel(0.),
  m_ch(0)
{
}
/*--------------------------------------------------------------------------------*/
// Fake Rate Destructor
/*--------------------------------------------------------------------------------*/
MeasureFakeRate2::~MeasureFakeRate2()
{
}
/*--------------------------------------------------------------------------------*/
// Begin -- Called at the beginning of the event loop
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::Begin(TTree* /*tree*/)
{
  if(m_dbg) cout << "MeasureFakeRate2::Begin" << endl;
  SusySelectionMatt::Begin(0);
  cout<<"Creating: "<<out<<endl;
  initHistos(m_fileName);
}
/*--------------------------------------------------------------------------------*/
// Terminate
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::Terminate()
{
  if(m_dbg) cout << "MeasureFakeRate2::Terminate" << endl;
  cout<<"Writing file "<<m_outFile<<endl;
  m_outFile->Write();
  cout<<"Closing file"<<endl;
  m_outFile->Close();
  dumpEventCounters();
}
/*--------------------------------------------------------------------------------*/
// Initialize the histograms
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::initHistos(string outName)
{
  // DG : here we will switch to a loop over m_controlRegions and m_signalRegions
  // const int CR_N(getNumControlRegions());
  if(CR_N >= kNmaxControlRegions)
    cout<<"MeasureFakeRate2::initHistos:"
        <<" "<<kNmaxControlRegions<<" histos reserved"
        <<" trying to book "<<CR_N<<endl
        <<" Exiting."
        <<endl;
  assert(CR_N<kNmaxControlRegions);
  cout<<"Creating file: "<<outName<<endl;
  m_outFile = new TFile((outName).c_str(),"recreate");
  m_outFile->cd();
  cout<<"File created: "<<m_outFile<<endl;

  // All of the rates will be stored in efficiency objects.
  // These are much more convenient for handling the errors.
  #define EFFVAR(eff,name,nbins,bins) \
    do{ eff = new EffObject(name,nbins,bins); } while(0)
  #define EFF(eff,name,nbins,xmin,xmax) \
    do{ eff = new EffObject(name,nbins,xmin,xmax); } while(0)
  #define LABEL(eff, bin, label) \
    do{ eff->SetXLabel(bin, label); }while(0)

  for(int il=0; il<LT_N; ++il){ // Loop over lepton type
    string lName = LTNames[il];
    // For the control regions we need the plots from both Data and MC
    for(int icr=0; icr<CR_N; ++icr){
      string cName = CRNames[icr];
      // Loop over the channels
      for(int ich=0; ich<Ch_N; ++ich){
        string chan = chanNames[ich];
        string base = lName + "_" + cName + "_" + chan + "_";
        EFFVAR(h_l_pt[il][icr][ich], (base+"l_pt"),nFakePtbins,FakePtbins);
        EFFVAR(h_l_pt_coarse[il][icr][ich], (base+"l_pt_coarse"),nCoarseFakePtbins,coarseFakePtbins);
        EFFVAR(h_l_pt_heavy[il][icr][ich], (base+"l_pt_heavy"),nFakePtbins,FakePtbins);
        EFFVAR(h_l_pt_others[il][icr][ich], (base+"l_pt_others"),nFakePtbins,FakePtbins);
        EFFVAR(h_l_eta[il][icr][ich], (base+"l_eta"),nEtabins,Etabins);
        EFFVAR(h_l_eta_coarse[il][icr][ich], (base+"l_eta_coarse"),nCoarseEtabins,CoarseEtabins);
        EFFVAR(h_metrel[il][icr][ich], (base+"metrel"),nMetbins, Metbins);
        EFFVAR(h_met[il][icr][ich], (base+"met"),nMetbins, Metbins);
        EFF(h_metrel_fine[il][icr][ich], (base+"metrel_fine"),40,0,200);
        EFF(h_met_fine[il][icr][ich], (base+"met_fine"),40,0,200);
        EFFVAR(h_metrel_coarse[il][icr][ich], (base+"metrel_coarse"),nCoarseMetbins, coarseMetbins);
        EFFVAR(h_met_coarse[il][icr][ich], (base+"met_coarse"),nCoarseMetbins, coarseMetbins);
        EFFVAR(h_njets[il][icr][ich], (base+"njets"),nJetbins, Jetbins);
        EFFVAR(h_nlightjets[il][icr][ich], (base+"nlightjets"),nJetbins, Jetbins);
        EFFVAR(h_nheavyjets[il][icr][ich], (base+"nheavyjets"),nJetbins, Jetbins);
        EFFVAR(h_nlightjetsNoB[il][icr][ich], (base+"nlightjetsNoB"),nJetbins, Jetbins);
        EFF(h_onebin[il][icr][ich], (base+"onebin"), 1, -0.5, 0.5);
        // d0sig
        EFF(h_heavy_d0sig[il][icr][ich], (base+"heavy_d0sig"), 50, 0, 5);
        EFF(h_light_d0sig[il][icr][ich], (base+"light_d0sig"), 50, 0, 5);
        EFF(h_conv_d0sig[il][icr][ich], (base+"conv_d0sig"), 50, 0, 5);
        EFF(h_l_type[il][icr][ich], (base+"l_type"), nType, Typemin, Typemax);
        EFF(h_l_origin[il][icr][ich], (base+"l_origin"), nOrigin, Originmin, Originmax);
        // Flavor counting
        EFF(h_flavor[il][icr][ich], (base+"flavor"), LS_N, -0.5, LS_N-0.5);
        for(int lbl=0; lbl<LS_N; ++lbl) LABEL(h_flavor[il][icr][ich], lbl+1, LSNames[lbl]);

        EFF(h_ht[il][icr][ich], (base+"ht"), 20, 0, 400);
        EFFVAR(h_ht_pt[il][icr][ich], (base+"ht_pt"), nHtbins, Htbins);
        EFF(h_ht_wMet[il][icr][ich], (base+"ht_wMet"), 20, 0, 400);
        EFFVAR(h_ht_pt_wMet[il][icr][ich], (base+"ht_pt_wMet"), nHtbins, Htbins);
        EFF(h_with_without_Etcone[il][icr][ich], (base+"with_without_Etcone"),3, -0.5, 2.5);
        }// end loop over channels
      // Saving CR plots to look at how cuts affect distribution
      string base = lName + "_" + cName + "_all_";
      EFF(h_met_cr[il][icr], (base+"met_cr"), 50, 0, 200);
      EFF(h_mt_tag_cr[il][icr], (base+"mt_tag_cr"), 50, 0, 200);
      EFF(h_mt_probe_cr[il][icr], (base+"mt_probe_cr"), 50, 0, 200);
      EFF(h_ht_cr[il][icr], (base+"ht_cr"), 50, 0, 400);
      EFF(h_mll_cr[il][icr], (base+"mll_cr"), 30, 0, 300);
    }// end loop over control regions
  }// end loop over lepton types
  #undef EFFVAR
  #undef EFF
  #undef LABEL
}
//----------------------------------------------------------
void printProgress(const Susy::Event *e, Long64_t counter)
{
  cout << "**** Processing entry " <<counter
       << " run "<<e->run
       << " event "<<e->event << " ****" << endl;
}
//----------------------------------------------------------
Bool_t MeasureFakeRate2::Process(Long64_t entry)
{
  if(m_dbg) cout << "MeasureFakeRate2::Process" << endl;
  GetEntry(entry);
  clearObjects();
  m_chainEntry++;
  increment(n_readin);
  if(m_dbg || m_chainEntry%50000==0) printProgress(nt.evt(), m_chainEntry);
  selectObjects(NtSys_NOM, false, TauID_medium);
  if( !selectEvent() ) return false;
  if(m_signalTaus.size()!=0)      /*increment(n_pass_CRWHSStauv  [m_ET]); else */ return false;
  if( m_baseLeptons.size() >= 2 ) increment(n_pass_atleast2Lep);
  if( m_baseLeptons.size() == 2 ) increment(n_pass_exactly2Lep);
  bool uniqueLepPair(m_baseLeptons.size() == 2);
  if(uniqueLepPair && !susy::passMllMin(m_baseLeptons, 20.0)) return false;
  // Loop over all regions (regardless of data or mc) and only fill relevant data quantitites.
  for(int cr = 0; cr<CR_N; ++cr){
    ControlRegion CR = (ControlRegion) cr;
    m_ch = Ch_all;
    m_probes.clear();
    m_tags.clear();
    bool passCR = false;
    const LeptonVector& leptons = m_baseLeptons;
    const JetVector& jets = m_signalJets2Lep;
    switch(CR) {
      // Data Driven CRs
    case CR_Real     : passCR = passRealCR  (leptons, jets,m_met, CR); break;
    case CR_SideLow  : passCR = passRealCR  (leptons, jets,m_met, CR); break;
    case CR_SideHigh : passCR = passRealCR  (leptons, jets,m_met, CR); break;
    case CR_HF       : passCR = passHFCR    (leptons, jets,m_met, CR); break;
    case CR_HF_high  : passCR = passHFCR    (leptons, jets,m_met, CR); break;
    case CR_LFZjet   : passCR = passLFZjetCR(leptons, jets,m_met    ); break;
    case CR_LFWjet   : passCR = passLFWjetCR(leptons, jets,m_met    ); break;
    case CR_Conv     : passCR = passConvCR  (leptons, jets,m_met    ); break;
    case CR_CFlip    : passCR = passCFCR    (leptons, jets,m_met    ); break;
      // MC Regions
    case CR_MCHeavy  : passCR = passMCReg(leptons,jets,m_met,CR); break;
    case CR_MCLight  : passCR = passMCReg(leptons,jets,m_met,CR); break;
    case CR_MCConv   : passCR = passMCReg(leptons,jets,m_met,CR); break;
    case CR_MCQCD    : passCR = passMCReg(leptons,jets,m_met,CR); break;
    case CR_MCALL    : passCR = passMCReg(leptons,jets,m_met,CR); break;
    case CR_MCReal   : passCR = passMCReg(leptons,jets,m_met,CR); break;
    case CR_MCNone   : passCR = passMCReg(leptons,jets,m_met,CR); break;
      // Remaining signal and control regions
    default          : passCR = passSignalRegion(leptons,jets,m_met,CR); break;
    } // end switch(CR)
    if( passCR ){
      for(uint ip=0; ip<m_probes.size(); ++ip)
        fillRatesHistos(m_probes.at(ip), jets, m_met, CR);
    }// end if Pass Control Region
  }// end loop over Control Regions
  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// Plotting Method
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::fillRatesHistos(const Lepton* lep, const JetVector& jets,
                                       const Met* met, ControlRegion CR)
{

  LeptonType lt(lep->isEle() ? LT_EL : LT_MU);
  bool pass(isSignalLepton(lep, m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC));

  #define FILL(eff, var) \
    do{                  eff[lt][CR][m_ch]->Fill(pass,m_evtWeight,var);   \
      if(m_ch != Ch_all) eff[lt][CR][Ch_all]->Fill(pass,m_evtWeight,var); \
    }while(0)

  FILL(h_l_pt, lep->Pt());
  FILL(h_l_pt_coarse, lep->Pt());
  FILL(h_l_eta, fabs(lep->Eta()));
  FILL(h_l_eta_coarse, fabs(lep->Eta()));
  if(nt.evt()->isMC){
    if( isHFLepton(lep) ) FILL(h_l_pt_heavy, lep->Pt());
    else                  FILL(h_l_pt_others, lep->Pt());
  }
  FILL(h_metrel, m_metRel);
  FILL(h_met, met->Et);
  FILL(h_metrel_coarse, m_metRel);
  FILL(h_met_coarse, met->Et);
  FILL(h_metrel_fine, m_metRel);
  FILL(h_met_fine, met->Et);
  // Number of jets
  FILL(h_njets, jets.size());
  FILL(h_nlightjets, numberOfCLJets(jets));
  FILL(h_nheavyjets, numberOfCBJets(jets));
  if( numberOfCBJets(jets) == 0 ) FILL(h_nlightjetsNoB, numberOfCLJets(jets));
  // Single bin
  FILL(h_onebin, 0);
  // If the event is MC, save the flavor
  if( nt.evt()->isMC ){
    LeptonSource ls = getLeptonSource(lep);
    if(ls == LS_Real){
      FILL(h_l_type, lep->mcType);
      FILL(h_l_origin, lep->mcOrigin);
    }
    FILL(h_flavor, ls);
    if(ls==LS_HF||ls==LS_LF) FILL(h_flavor, LS_QCD);
    if(ls == LS_HF)   FILL(h_heavy_d0sig, lep->d0Sig(true));
    if(ls == LS_LF)   FILL(h_light_d0sig, lep->d0Sig(true));
    if(ls == LS_Conv) FILL(h_conv_d0sig, lep->d0Sig(true));
  }
  else{
    float d0sig = lep->d0Sig(true);
    FILL(h_heavy_d0sig, d0sig);
    FILL(h_light_d0sig, d0sig);
    FILL(h_conv_d0sig, d0sig);
  }

  // Ht Plots
  float ht = lep->Pt();
  for(uint ij=0; ij<jets.size(); ++ij) ht += jets.at(ij)->Pt();
  FILL(h_ht, ht);
  FILL(h_ht_pt, ht/lep->Pt());
  FILL(h_ht_wMet, ht+met->Et);
  FILL(h_ht_pt_wMet, (ht+met->Et)/lep->Pt());

  bool with = isSignalWithEtcone(lep);
  bool without = isSignalWithoutEtcone(lep);
  FILL(h_with_without_Etcone, 0);
  if(with)    FILL(h_with_without_Etcone, 1);
  if(without) FILL(h_with_without_Etcone, 2);
  #undef FILL

}
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::fillCrHistos(const Lepton* tag, const Lepton* probe, const JetVector& jets, const Met* met,
                                    ControlRegion CR, CRPLOT CRP)
{
  LeptonType lt = probe->isEle() ? LT_EL : LT_MU;
  bool pass(isSignalLepton(probe, m_baseElectrons, m_baseMuons, nt.evt()->nVtx, nt.evt()->isMC));
  #define FILL(eff, var)							\
    do{ eff[lt][CR]->Fill(pass,m_evtWeight,var); } while(0)

  if(CRP == CRP_mll)       FILL(h_mll_cr, (*tag+*probe).M());
  if(CRP == CRP_met)       FILL(h_met_cr, met->Et);
  if(CRP == CRP_mt_tag)    FILL(h_mt_tag_cr, Mt(tag, met));
  if(CRP == CRP_mt_probe)  FILL(h_mt_probe_cr, Mt(probe, met));
  if(CRP == CRP_ht){
    float ht = tag->Pt() + probe->Pt();
    for(uint i=0; i<jets.size(); ++i) ht += jets.at(i)->Pt();
    FILL(h_ht_cr, ht);
  }
  #undef FILL
}

/*--------------------------------------------------------------------------------*/
// Control Regions
/*--------------------------------------------------------------------------------*/

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// MC Control region is defined by metrel
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passMCReg(const LeptonVector &leptons,
				 const JetVector &jets,
				 const Met* met,
				 ControlRegion CR)
{
  // Used to measure fake rate in MC in a region
  // that is close to our signal region
  // * 40 < metrel < 100
  // * exactly two leptons
  // * probes are the fake leptons of specific type
  if( !nt.evt()->isMC )     return false;
  if( leptons.size() != 2 ) return false;
  m_ch = getChan(leptons);
  m_metRel = getMetRel(met,leptons,jets);
  for(uint il=0; il<leptons.size(); ++il){
    Lepton* l=leptons[il];
    if(CR==CR_MCHeavy && isHFLepton(l))  m_probes.push_back( l ); // Heavy
    if(CR==CR_MCLight && isLFLepton(l))  m_probes.push_back( l ); // Light
    if(CR==CR_MCConv && isConvLepton(l)) m_probes.push_back( l ); // Conversion
    if(CR==CR_MCNone && isFakeLepton(l)) m_probes.push_back( l ); // None CR
    if(CR==CR_MCQCD && isQCDLepton(l))   m_probes.push_back( l ); // QCD
    if(CR==CR_MCALL)                     m_probes.push_back( l ); // All fakes <-- Will include reals
    if(CR==CR_MCReal && isRealLepton(l)) m_probes.push_back( l ); // Real
  }
  m_evtWeight = getEvtWeight(leptons);
  return true;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// The signal region checks
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passSignalRegion(const LeptonVector &leptons,
                                        const JetVector &jets,
                                        const Met* met,
                                        ControlRegion CR)
{
  if( leptons.size() != 2 ) return false;
  bool passSR = false;
  switch(CR) {
  case CR_SRmT2a       : passSR = passSRmT2a    (leptons,jets,met,false,true); break;
  case CR_SRmT2b       : passSR = passSRmT2b    (leptons,jets,met,false,true); break;
  case CR_SRmT2c       : passSR = passSRmT2c    (leptons,jets,met,false,true); break;
  case CR_SRWWa        : passSR = passSRWWa     (leptons,jets,met,false,true); break;
  case CR_SRWWb        : passSR = passSRWWb     (leptons,jets,met,false,true); break;
  case CR_SRWWc        : passSR = passSRWWc     (leptons,jets,met,false,true); break;
  case CR_SRZjets      : passSR = passSRZjets   (leptons,jets,met,false,true); break;
  case CR_VRSS         : passSR = passVRSS      (leptons,jets,met);            break;
  case CR_SSInc        : passSR = sameSign      (leptons);                     break;
  case CR_CRWWMet      : passSR = passCRWWMet   (leptons,jets,met,false,true); break;
  case CR_CRWWmT2      : passSR = passCRWWmT2   (leptons,jets,met,false,true); break;
  case CR_CRTopMet     : passSR = passCRTopMet  (leptons,jets,met,false,true); break;
  case CR_CRTopmT2     : passSR = passCRTopmT2  (leptons,jets,met,false,true); break;
  case CR_CRZVMet      : passSR = passCRZVMet   (leptons,jets,met,false,true); break;
  case CR_CRZVmT2_90   : passSR = passCRZVmT2a  (leptons,jets,met,false,true); break;
  case CR_CRZVmT2_120  : passSR = passCRZVmT2b  (leptons,jets,met,false,true); break;
  case CR_CRZVmT2_150  : passSR = passCRZVmT2c  (leptons,jets,met,false,true); break;
  case CR_CRZVmT2_100  : passSR = passCRZVmT2d  (leptons,jets,met,false,true); break;
  case CR_CRTopZjets   : passSR = passCRTopZjets(leptons,jets,met,false,true); break;
  case CR_CRZXZjets    : passSR = passCRZXZjets (leptons,jets,met,false,true); break;
  case CR_PremT2       : passSR = passCRPremT2  (leptons,jets,met);            break;
  case CR_SRWHSS       : passSR = passWhSS      (leptons,jets,met);            break;
  default: cout<<"invalid ControlRegion "<<CR<<endl;
  }
  for(uint i=0; i<leptons.size(); ++i) m_probes.push_back( leptons[i]);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);
  m_ch = getChan(leptons);
  return passSR;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Real Control Region: Z
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passRealCR(const LeptonVector &leptons,
                                  const JetVector &jets,
                                  const Met* met,
                                  ControlRegion CR)
{
  // Real CRA:
  // * Require Exactly two baseline leptons SF
  // * Require single lepton triggers
  // * Require fall in Z mass window
  // * Require at least one tight lepton for tag

  // Same flavor dilepton
  if( leptons.size() != 2 )  return false;
  if( !sameFlavor(leptons) ) return false;
  // In Z window or side bands
  float mll = Mll(leptons[0],leptons[1]);
  if(CR == CR_Real)         { if( fabs(mll-91.2) > 10 ) return false; }
  else if(CR == CR_SideLow) { if( !(61 < mll && mll < 71) )   return false; }
  else if(CR == CR_SideHigh){ if( !(111 < mll && mll < 121) ) return false; }
  // At least one tight lepton
  int nVtx   = nt.evt()->nVtx;
  bool isMC  = nt.evt()->isMC;
  bool l0sig(isSignalLepton(leptons[0], m_baseElectrons, m_baseMuons, nVtx, isMC));
  bool l1sig(isSignalLepton(leptons[1], m_baseElectrons, m_baseMuons, nVtx, isMC));
  if( !l0sig && !l1sig ) return false;
  // Pass trigger
  uint l0flag(leptons[0]->trigFlags), l1flag(leptons[1]->trigFlags);
  bool l0trig(leptons[0]->isEle() ? l0flag & TRIG_e24vhi_medium1 : l0flag & TRIG_mu24i_tight);
  bool l1trig(leptons[1]->isEle() ? l1flag & TRIG_e24vhi_medium1 : l1flag & TRIG_mu24i_tight);
  if( l0sig && l0trig && leptons[0]->Pt() > 25 ) m_probes.push_back( leptons[1] );
  if( l1sig && l1trig && leptons[1]->Pt() > 25 ) m_probes.push_back( leptons[0] );
  if( m_probes.size() == 0 ) return false;
  m_metRel = getMetRel(met, leptons, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);
  return true;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Heavy Flavor Control region
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passHFCR(const LeptonVector &leptons,
                                const JetVector &jets,
                                const Met* met,
                                ControlRegion CR)
{
  // This is a heavy flavor tag and probe trying
  // to select b-bbar events by looking at loose
  // muon tag in the event.
  // * Require Muon Stream Only
  // * Exactly one muon overlapping exactly one b-jet
  // * Tag muon Pt > 20 and match to dilepton trigger
  //    -- EF_mu18_tight_e7_medium (elec probe)
  //    -- EF_mu18_tight_mu8_EFFS  (muon probe)
  // * Probe lepton baseline after overlap
  // * Mt(probe,met) < 40 GeV
  // * met < 40 GeV
  // * Veto dimuon events 10 GeV near Z
  // * m(ll) > 40 GeV (probe is muon)

  // Check pre muons and jets
  MuonVector preMuons = getPreMuons(&nt, NtSys_NOM);
  if( preMuons.size() == 0 ) return false;
  if( jets.size() == 0 )     return false;

  Lepton* tag = NULL; // Find tag
  int nTags   = 0;
  int nBjets  = 0;
  for(uint ij = 0; ij<jets.size(); ++ij){
    if( !isCentralBJet(jets.at(ij)) ) continue;
    nBjets++;
    for(uint imu=0; imu<preMuons.size(); ++imu){
      if(jets.at(ij)->DeltaR( *preMuons.at(imu) ) < 0.4){
        tag = (Lepton*) preMuons.at(imu);
        nTags++;
      }// end if overlap
    }// end loop over muons
  }// end loop over jets
  if( nTags != 1 )  return false;
  if( nBjets != 1 ) return false;
  Lepton* probe = NULL; // Get the probe
  uint nProbes  = 0;
  for(uint il=0; il<leptons.size(); ++il){
    if     (leptons.at(il)->isEle()){ probe = leptons.at(il); nProbes++; }
    else if( tag != leptons.at(il) ){ probe = leptons.at(il); nProbes++; }
  }// end loop over baseline leptons
  if( nProbes != 1 ) return false;
  // Check mll
  fillCrHistos(tag,probe,jets,met,CR_HF,CRP_mll);
  bool isMM = tag->isMu() && probe->isMu();
  if(isMM){
    float mll(Mll(tag,probe)), mZ(91.2);
    if( fabs(mll-mZ) < 10 ) return false;
    if( mll < 40 )          return false;
  }
  // Check the trigger
  if( tag->Pt() < 20 ) return false;
  uint tagFlag = tag->trigFlags;
  if( !isMM && !((tagFlag & TRIG_mu18_tight) && (tagFlag & TRIG_mu18_tight_e7_medium1)) )
    return false;
  if( isMM && !((tagFlag & TRIG_mu18_tight) && (tagFlag & TRIG_mu18_tight_mu8_EFFS)) )
    return false;

  // Met and Mt Cut
  fillCrHistos(tag,probe,jets,met,CR_HF,CRP_ht);
  fillCrHistos(tag,probe,jets,met,CR_HF,CRP_met);
  if( met->Et > 40 )       return false;
  fillCrHistos(tag,probe,jets,met,CR_HF,CRP_mt_tag);
  fillCrHistos(tag,probe,jets,met,CR_HF,CRP_mt_probe);
  if( CR == CR_HF )      if( Mt(probe,met) > 40  ) return false;
  if( CR == CR_HF_high ) if( Mt(probe,met) > 100 ) return false;

  // We have survived so save objects and return
  LeptonVector temp; temp.push_back(tag); temp.push_back(probe);
  m_metRel = getMetRel(met, temp, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(temp, true);
  m_probes.push_back( probe );
  m_tags.push_back( tag );
  return true;
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Light Flavor Control Region -- Zjet
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passLFZjetCR(const LeptonVector &leptons,
				    const JetVector &jets,
				    const Met* met)
{
  // Light flavor tag and probe using Z+jet events
  // * Same-flavor OS signal leptons 10 GeV of Z peak
  // * Probe lepton of different flavor
  // * M(lll) not within 10 GeV of Z peak
  // * met < 40 GeV
  // * Mt(probe,met) < 50 GeV
  // * Reject events with b-jet
  // * Pass LF Trig function

  // Require exactly 3 leptons
  if( leptons.size() != 3 )       return false;

  // Reject events with b-jet
  if( numberOfCBJets(jets) != 0 ) return false;

  // Use Signal info since it's cleaner
  LeptonVector sigLeps;
  if(m_signalMuons.size() == 2 && m_baseElectrons.size() == 1){
    if(m_signalMuons[0]->q * m_signalMuons[1]->q < 0){
      sigLeps.push_back( (Lepton*) m_signalMuons[0]);
      sigLeps.push_back( (Lepton*) m_signalMuons[1]);
      m_probes.push_back( (Lepton*) m_baseElectrons[0]);
    }
  }
  if(m_baseMuons.size() == 1 && m_signalElectrons.size() == 2){
    if(m_signalElectrons[0]->q * m_signalElectrons[1]->q < 0){
      sigLeps.push_back( (Lepton*) m_signalElectrons[0]);
      sigLeps.push_back( (Lepton*) m_signalElectrons[1]);
      m_probes.push_back( (Lepton*) m_baseMuons[0]);
    }
  }
  if(sigLeps.size() != 2)    return false;
  if( !passLFTrig(sigLeps) ) return false;
  // Mass rejections
  float mZ   = 91.2;
  float mll  = Mll(sigLeps[0],sigLeps[1]);
  float mlll = Mlll(sigLeps[0],sigLeps[1],m_probes[0]);
  if( fabs(mll-mZ) > 10 )  return false;
  if( fabs(mlll-mZ) < 10 ) return false;
  // Mt Cut
  if( Mt(m_probes[0],met) > 50 ) return false;
  // Update and save objects
  m_metRel = getMetRel(met, leptons, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(sigLeps, true) * m_probes[0]->effSF;
  return true;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Light Flavor Control Region -- Wjet
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passLFWjetCR(const LeptonVector &leptons,
				    const JetVector &jets,
				    const Met* met)
{

  // Light Flavor control region for W+jet
  // * Two same-sign same-flavor leptons
  // * Tag is leading probe is subleading
  // * Tag must be signal and match to tight trigger
  // * Veto mll region 70-100 get rid of Zs
  // * Mt(tag,met) > 40
  // * Met/HT > 0.4
  // * Met > 20 GeV
  // * No b jets

  // Select two baseline leptons of SF SS
  // Reject events with bJet
  if( leptons.size() != 2 )     return false;
  if( oppositeFlavor(leptons) ) return false;
  if( oppositeSign(leptons) )   return false;
  // Define tag as leading lepton. This is very nearly true for W+jet
  Lepton* tag = leptons[0];
  Lepton* probe = leptons[1];
  // Require tag to be signal lepton
  if( !isSignalLepton(tag, m_baseElectrons, m_baseMuons, nt.evt()->nVtx, nt.evt()->isMC) )
    return false;
  // Make sure lepton passes single trigger
  if( tag->isEle() ){
    if( !(tag->Pt() > 25 && tag->matchTrig(TRIG_e24vhi_medium1)) ) return false;
  }
  else if(tag->isMu()){
    if( !(tag->Pt() > 25 && tag->matchTrig(TRIG_mu24i_tight)) ) return false;
  }
  // First cut: Reject Z
  float mll = Mll(tag,probe);
  fillCrHistos(tag,probe,jets,met,CR_LFWjet,CRP_mll);
  if( 70 < mll && mll < 100 ) return false;
  // Second Cut: Mt to reduce Z and enhace W+jet
  fillCrHistos(tag,probe,jets,met,CR_LFWjet,CRP_mt_tag);
  fillCrHistos(tag,probe,jets,met,CR_LFWjet,CRP_mt_probe);
  if( Mt(tag,met) < 40 ) return false;
  // Thid cut: MET/HT shown to separate W+jet and all other backgrounds
  float Ht = tag->Pt() + probe->Pt();
  for(uint i=0; i<jets.size(); ++i) Ht += jets.at(i)->Pt();
  fillCrHistos(tag,probe,jets,met,CR_LFWjet,CRP_ht);
  if( met->Et/Ht < 0.4 ) return false;
  // Fourth Cut: Met > 20 to get rid of region with bad agreement..
  fillCrHistos(tag,probe,jets,met,CR_LFWjet,CRP_met);
  if( met-> Et < 20 ) return false;
  // We made it!
  m_metRel = getMetRel(met,leptons,jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons, true);
  m_probes.push_back(probe);
  return true;
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Conversion Control Region -- Conv
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passConvCR(const LeptonVector &leptons,
                                  const JetVector &jets,
                                  const Met* met)
{
  // Conversion control region to get electron conversion rate
  // * Data only from Muon Stream
  // * Exactly two signal muons of opposite sign
  // * Exactly one baseline electron
  // * M(lll) in enlarged Z window 80-100
  // * Met < 50 GeV
  // * M(mu,mu) > 20
  // * Mt(elec, met) < 40
  ElectronVector preElecs; MuonVector preMuons;
  if( m_signalMuons.size() != 2 )  return false;
  if( m_baseElectrons.size() != 1) return false;
  if( m_baseLeptons.size() != 3)   return false;
  preMuons.push_back( m_signalMuons[0] );
  preMuons.push_back( m_signalMuons[1] );
  preElecs.push_back( m_baseElectrons[0] );
  float mlll = (*preMuons[0] + *preMuons[1] + *preElecs[0]).M();
  if( !(80 < mlll && mlll < 100) )           return false;
  if( preMuons[0]->q * preMuons[1]->q > 0 )  return false; // Opposite sign
  LeptonVector tempL; ElectronVector tempE;
  buildLeptons(tempL, tempE, preMuons);
  if( !passTrigger(tempL) )                  return false; // Pt and trigger requirement
  if( met->Et > 50 )                         return false; // Met < 50
  if( (*preMuons[0]+*preMuons[1]).M() < 20 ) return false; // M(mu,mu) > 20
  if( Mt((Lepton*) preElecs[0], met) > 40)   return false; // Mt(elec, met) < 40
  m_probes.push_back( (Lepton*) preElecs[0] );
  m_metRel = getMetRel(met, leptons, jets); // Not useful quantity for this region.
  float eff = nt.evt()->isMC ? m_probes[0]->effSF : 1.0;
  m_evtWeight = getEvtWeight(tempL,true,true) * eff;
  return true;
}

/*--------------------------------------------------------------------------------*/
// Pass Charge flip control region
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::passCFCR(const LeptonVector &leptons,
                                const JetVector& jets,
                                const Met* met)
{
  // Real CRA:
  // * Require Exactly two baseline leptons SF
  // * Require single lepton triggers
  // * Require fall in Z mass window
  // * Require at least one tight lepton for tag
  // * Require to be same-sign

  // Same flavor dilepton
  if( leptons.size() != 2 )  return false;
  if( !sameFlavor(leptons) ) return false;
  if( !sameSign(leptons) )   return false;
  // In Z window or side bands
  float mll = Mll(leptons[0],leptons[1]);
  if( !(70 < mll && mll < 100) ) return false;
  // At least one tight lepton
  int nVtx   = nt.evt()->nVtx;
  bool isMC  = nt.evt()->isMC;
  bool l0sig(isSignalLepton(leptons[0], m_baseElectrons, m_baseMuons, nVtx, isMC));
  bool l1sig(isSignalLepton(leptons[1], m_baseElectrons, m_baseMuons, nVtx, isMC));
  if( !l0sig && !l1sig ) return false;
  // Pass trigger
  uint l0flag = leptons[0]->trigFlags;
  bool l0trig = leptons[0]->isEle() ? l0flag & TRIG_e24vhi_medium1 : l0flag & TRIG_mu24i_tight;
  uint l1flag = leptons[1]->trigFlags;
  bool l1trig = leptons[1]->isEle() ? l1flag & TRIG_e24vhi_medium1 : l1flag & TRIG_mu24i_tight;
  if( l0sig && l0trig && leptons[0]->Pt() > 25 ) m_probes.push_back( leptons[1] );
  if( l1sig && l1trig && leptons[1]->Pt() > 25 ) m_probes.push_back( leptons[0] );
  if( m_probes.size() == 0 ) return false;
  m_metRel = getMetRel(met, leptons, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);
  return true;
}

/*--------------------------------------------------------------------------------*/
// Light Flavor Trigger Requirements
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::passLFTrig(const LeptonVector& leps)
{
  // Try LF control region with the same trigger scheme
  // that 3-lepton uses
  if( leps.size() != 2 ) return false;
  uint ie=0;
  uint im=0;
  for(uint i=0; i<leps.size(); ++i){ if(leps[0]->isEle()) ie++; else im++; }
  bool passMuons = false;
  uint n1M = 0;
  uint nSym2M = 0;
  uint nAsym2M = 0;
  uint nAsym2M_m18 = 0;
  // Trigger matching for Muon Stream
  for(uint i=0; i<leps.size(); i++){
    if(leps.at(i)->isEle()) continue;
    const Muon* lep = (Muon*) leps[i];
    float pt = lep->Pt();
    // Single muon trigger
    if(pt>25 && lep->matchTrig(TRIG_mu24i_tight)) n1M++;
    // 2m symmetric trigger
    if(pt>14 && lep->matchTrig(TRIG_2mu13)) nSym2M++;
    // 2m asymmetric trigger
    if(lep->matchTrig(TRIG_mu18_tight_mu8_EFFS)){
      nAsym2M++;
      if(pt>18 &&  lep->matchTrig(TRIG_mu18_tight)) nAsym2M_m18++;
    }
  }
  // Logical OR or all Muon Triggers
  if(n1M || nSym2M || nAsym2M || nAsym2M_m18) passMuons = true;
  bool passEgamma = false;
  uint n1E = 0;
  uint nSym2E = 0;
  uint nAsym2E = 0;
  uint nAsym2E_e24 =0;
  // Trigger matching for Electron Stream
  for(uint i=0; i<leps.size(); i++){
    if(!leps.at(i)->isEle()) continue;
    const Electron* lep = (Electron*) leps[i];
    float pt = lep->Pt();
    // Single electron trigger
    if(pt>25 && lep->matchTrig(TRIG_e24vhi_medium1)) n1E++;
    // 2e symmetric trigger
    if(pt>14 && lep->matchTrig(TRIG_2e12Tvh_loose1)) nSym2E++;
    // 2e asymmetric trigger
    if(lep->matchTrig(TRIG_e24vh_medium1_e7_medium1)){
      nAsym2E++;
      if(pt>25 && lep->matchTrig(TRIG_e24vh_medium1)) nAsym2E_e24++;
    }
  }
  // Logical OR or all Electron triggers
  if(n1E || nSym2E || nAsym2E || nAsym2E_e24) passEgamma = true;
  if( !nt.evt()->isMC ){
    DataStream stream = nt.evt()->stream;
    if( im == 2 && stream == Stream_Muons) return passMuons;
    if( ie == 2 && stream == Stream_Egamma) return passEgamma;
    return false;
  }
  return (passMuons && im==2) || (passEgamma && ie ==2);
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithEtcone(const Lepton* lep)
{
  bool isSignal = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  // If Electron, just return
  if( lep->isEle() ) return isSignal;
  // Otherwise muon, add etcone cut:
  return isSignal && ((Muon*) lep)->etcone30/lep->Pt() < 0.18;
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutEtcone(const Lepton* lep)
{
  // If Electron, just return
  if( lep->isMu() )
    return isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  // Otherwise, turn off for electrons
  m_doElEtconeCut = false;
  bool isSignal = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doElEtconeCut = true;
  return isSignal;
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutPtcone(const Lepton* lep)
{
  m_doPtconeCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doPtconeCut = true;
  return pass;
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutd0Sig(const Lepton* lep)
{
  m_doIPCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doIPCut = true;
  float cutVal = lep->isEle() ? ELECTRON_Z0_SINTHETA_CUT : MUON_Z0_SINTHETA_CUT;
  bool passz0  = fabs( lep->z0SinTheta(true) ) < cutVal;
  return pass && passz0;
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutz0Sig(const Lepton* lep)
{
  m_doIPCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doIPCut = true;
  float cutVal = lep->isEle() ? ELECTRON_D0SIG_CUT : MUON_D0SIG_CUT;
  bool passd0  = fabs( lep->d0Sig(true) ) < cutVal;
  return pass && passd0;
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutIP(const Lepton* lep)
{
  m_doIPCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doIPCut = true;
  return pass;
}
//----------------------------------------------------------
DiLepPair MeasureFakeRate2::getDilepPair(const Lepton* tag, const Lepton* probe){
  if( !tag ) return DL_LL; // no tag
  int nvtx  = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  bool lead = tag->Pt() >= probe->Pt() ?
    isSignalLepton(tag,m_baseElectrons,m_baseMuons,nvtx,isMC) :
    isSignalLepton(probe,m_baseElectrons,m_baseMuons,nvtx,isMC);
  bool sublead = tag->Pt() < probe->Pt() ?
    isSignalLepton(tag,m_baseElectrons,m_baseMuons,nvtx,isMC) :
    isSignalLepton(probe,m_baseElectrons,m_baseMuons,nvtx,isMC);
  if( lead && sublead )  return DL_TT;
  if( lead && !sublead ) return DL_TL;
  if( !lead && sublead ) return DL_LT;
  return DL_LL;
}
//----------------------------------------------------------
LeptonSource MeasureFakeRate2::getLeptonSource(const Lepton* l)
{
  if( isRealLepton(l) ) return LS_Real;
  if( isHFLepton(l) )   return LS_HF;
  if( isLFLepton(l) )   return LS_LF;
  if( isConvLepton(l) ) return LS_Conv;
  return LS_Unk;
}
