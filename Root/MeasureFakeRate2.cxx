////////////////////////////////////////////////////
// Will house the control regions used to measure //
// the fake rates used in 2-lep analysis.         //
////////////////////////////////////////////////////

#include "SusyTest0/MeasureFakeRate2.h"
#include "SusyTest0/SelectionRegions.h"
#include "SusyTest0/criteria.h"
#include "SusyTest0/FakeBinnings.h"

using namespace susywh; // pull in SelectionRegions
const int controlRegions[] = {
  kCR_Real, kCR_SideLow, kCR_SideHigh, kCR_HF, kCR_HF_high,
  kCR_Conv,
  kCR_MCConv, kCR_MCQCD,
  kCR_MCReal
};
const size_t nControlRegions = sizeof(controlRegions)/sizeof(controlRegions[0]);
const int signalRegions[] = {
  kCR_SSInc,
  kCR_SRWHSS,
  kPR_CR8lpt   ,
  kPR_CR8ee    ,
  kPR_CR8mm    ,
  kPR_CR8mmMtww,
  kPR_CR8mmHt
};
const size_t nSignalRegions = sizeof(signalRegions)/sizeof(signalRegions[0]);

/*--------------------------------------------------------------------------------*/
// Fake Rate Constructor
/*--------------------------------------------------------------------------------*/
MeasureFakeRate2::MeasureFakeRate2() :
  CR_N(nControlRegions + nSignalRegions),
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
    using namespace susy::fake; // pull in binning(s)
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
  // All the rates will be stored in efficiency objects.
  #define EFFVAR(eff,name,nbins,bins) do{ eff = new EffObject(name,nbins,bins); } while(0)
  #define EFF(eff,name,nbins,xmin,xmax) do{ eff = new EffObject(name,nbins,xmin,xmax); } while(0)
  for(int il=0; il<LT_N; ++il){ // for lName in [elec, muon]
    string lepton = LTNames[il];
    for(int icr=0; icr<CR_N; ++icr){ // see CRNames in SusyAnaDefsMatt.h for complete list
      string region = CRNames[icr];
      for(int ich=0; ich<Ch_N; ++ich){ // for chan in [all, ee, mm, em]
        string channel = chanNames[ich];
        string bn(lepton+"_"+region+"_"+channel+"_");
        h_l_pt         [il][icr][ich] = new EffObject(bn+"l_pt",         nFakePtbins,       FakePtbins);
        h_l_pt_coarse  [il][icr][ich] = new EffObject(bn+"l_pt_coarse",  nCoarseFakePtbins, coarseFakePtbins);
        h_l_eta        [il][icr][ich] = new EffObject(bn+"l_eta",        nEtabins,          Etabins);
        h_l_eta_coarse [il][icr][ich] = new EffObject(bn+"l_eta_coarse", nCoarseEtabins,    CoarseEtabins);
        h_metrel       [il][icr][ich] = new EffObject(bn+"metrel",       nMetbins,          Metbins);
        h_met          [il][icr][ich] = new EffObject(bn+"met",          nMetbins,          Metbins);
        h_njets        [il][icr][ich] = new EffObject(bn+"njets",        nJetbins,          Jetbins);
        h_onebin       [il][icr][ich] = new EffObject(bn+"onebin",       1,    -0.5, 0.5);
        h_flavor       [il][icr][ich] = new EffObject(bn+"flavor",       LS_N, -0.5, LS_N-0.5);
        for(int lbl=0; lbl<LS_N; ++lbl) h_flavor[il][icr][ich]->SetXLabel(lbl+1, LSNames[lbl]);
      } // end for(ich)
    } // end for(icr)
  } // end for(il)
  #undef EFFVAR
  #undef EFF
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
    ControlRegion CR = static_cast<ControlRegion>(cr);
    m_ch = Ch_all;
    m_probes.clear();
    m_tags.clear();
    bool passCR = false;
    const LeptonVector& leptons = m_baseLeptons;
    const JetVector& jets = m_signalJets2Lep;
    switch(CR) {
    case CR_Real     : passCR = passRealCR  (leptons, jets, m_met, CR); break;
    case CR_SideLow  : passCR = passRealCR  (leptons, jets, m_met, CR); break;
    case CR_SideHigh : passCR = passRealCR  (leptons, jets, m_met, CR); break;
    case CR_HF       : passCR = passHFCR    (leptons, jets, m_met, CR); break;
    case CR_HF_high  : passCR = passHFCR    (leptons, jets, m_met, CR); break;
    case CR_Conv     : passCR = passConvCR  (leptons, jets, m_met    ); break;
    case CR_MCConv   : passCR = passMCReg   (leptons, jets, m_met, CR); break;
    case CR_MCQCD    : passCR = passMCReg   (leptons, jets, m_met, CR); break;
    case CR_MCReal   : passCR = passMCReg   (leptons, jets, m_met, CR); break;
        // Remaining signal and control regions
    default          : passCR = passSignalRegion(leptons,jets,m_met,CR); break;
    } // end switch(CR)
    if( passCR ){
      for(size_t ip=0; ip<m_probes.size(); ++ip) fillRatesHistos(m_probes.at(ip), jets, m_met, CR);
    } // if(passCR)
  } // for(cr)
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

  FILL(h_l_pt,         lep->Pt());
  FILL(h_l_pt_coarse,  lep->Pt());
  FILL(h_l_eta,        fabs(lep->Eta()));
  FILL(h_l_eta_coarse, fabs(lep->Eta()));
  FILL(h_metrel,        m_metRel);
  FILL(h_met,           met->Et);
  FILL(h_njets,      jets.size());
  FILL(h_onebin, 0);  
  if( nt.evt()->isMC ){ // If the event is MC, save the flavor
    LeptonSource ls(getLeptonSource(lep));
    FILL(h_flavor, ls);
    if(ls==LS_HF||ls==LS_LF) FILL(h_flavor,      LS_QCD);
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
    if(CR==CR_MCConv && isConvLepton(l)) m_probes.push_back( l ); // Conversion
    if(CR==CR_MCQCD && isQCDLepton(l))   m_probes.push_back( l ); // QCD
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
  bool computeWhssPass(CR==CR_SRWHSS || CR==CR_CR8lpt || CR==CR_CR8ee || CR==CR_CR8mm
                       ||CR==CR_CR8mmMtww||CR==CR_CR8mmHt);
  SsPassFlags whssFlags(computeWhssPass ? passWhSS(leptons,jets,met) : SsPassFlags());
  bool isEe(Ch_ee==getChan(leptons)), isMm(Ch_mm==getChan(leptons));
  switch(CR) {
  case CR_SSInc        : passSR = sameSign(leptons);                           break;
  case CR_SRWHSS       : passSR =  whssFlags.metrel;                           break;
  case CR_CR8lpt       : passSR =  whssFlags.lepPt;                            break;
  case CR_CR8ee        : passSR = (whssFlags.lepPt   && isEe);                 break;
  case CR_CR8mm        : passSR = (whssFlags.lepPt   && isMm);                 break;
  case CR_CR8mmMtww    : passSR = (whssFlags.mtllmet && isMm);                 break;
  case CR_CR8mmHt      : passSR = (whssFlags.ht      && isMm);                 break;
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
  if( met->Et > 40 )       return false;
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
