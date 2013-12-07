#include "SusyTest0/MeasureFakeRate2.h"
#include "SusyTest0/criteria.h"
#include "SusyTest0/FakeBinnings.h"
#include "SusyTest0/SusySelection.h" // passEwkSs*

#include <algorithm> // transform
using namespace susy::fake;
using namespace susy::wh;
namespace sf = susy::fake;

const sf::Region controlRegions[] = {
    sf::CR_Real, sf::CR_SideLow, sf::CR_SideHigh, sf::CR_HF, sf::CR_HF_high,
    sf::CR_Conv,
    sf::CR_MCConv, sf::CR_MCQCD,
    sf::CR_MCReal
};
const size_t nControlRegions = sizeof(controlRegions)/sizeof(controlRegions[0]);
const sf::Region signalRegions[] = {
  sf::CR_SSInc,
  sf::CR_SRWHSS,
  sf::CR_CR8lpt,
  sf::CR_CR8ee,
  sf::CR_CR8mm,
  sf::CR_CR8mmMtww,
  sf::CR_CR8mmHt,
  sf::CR_SsEwk,
  sf::CR_SsEwkLoose
};
const size_t nSignalRegions = sizeof(signalRegions)/sizeof(signalRegions[0]);
const MeasureFakeRate2::LeptonType leptonTypes[] = {MeasureFakeRate2::kElectron, MeasureFakeRate2::kMuon};
const size_t nLeptonTypes(sizeof(leptonTypes) / sizeof(leptonTypes[0]));

/*--------------------------------------------------------------------------------*/
// Fake Rate Constructor
/*--------------------------------------------------------------------------------*/
MeasureFakeRate2::MeasureFakeRate2() :
  m_controlRegions(controlRegions, controlRegions + nControlRegions),
  m_signalRegions (signalRegions,  signalRegions  + nSignalRegions),
  m_leptonTypes(leptonTypes, leptonTypes + nLeptonTypes),
  m_fileName("measurefakerate2.root"),
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
  size_t nRegions(allRegions().size()), nLeptonTypes(m_leptonTypes.size());
  if(nRegions > kNmaxControlRegions)
    cout<<" Trying book histos for "<<nRegions<<" regions "<<" >= "<<kNmaxControlRegions<<endl<<" Exiting."<<endl;
  if(nLeptonTypes > kNmaxLeptonTypes)
    cout<<" Trying book histos for "<<nLeptonTypes<<" lepton types >= "<<kNmaxLeptonTypes<<endl<<" Exiting."<<endl;
  assert(nRegions <= kNmaxControlRegions);
  assert(nLeptonTypes <= kNmaxLeptonTypes);
  cout<<"Creating file: "<<outName<<endl;
  m_outFile = new TFile((outName).c_str(),"recreate");
  m_outFile->cd();
  cout<<"File created: "<<m_outFile<<endl;
  for(size_t il=0; il<m_leptonTypes.size(); ++il){
    string lepton(LeptonType2str(m_leptonTypes[il]));
    vector<sf::Region> regions(allRegions());
    for(size_t icr=0; icr<regions.size(); ++icr){
      string region(sf::region2str(regions[icr]));
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
  vector<sf::Region> regions(allRegions());
  for(size_t cr = 0; cr<regions.size(); ++cr){
    sf::Region CR = regions[cr];
    m_ch = Ch_all;
    m_probes.clear();
    m_tags.clear();
    bool passCR = false;
    const LeptonVector& leptons = m_baseLeptons;
    const JetVector& jets = m_signalJets2Lep;
    switch(CR) {
    case sf::CR_Real     : passCR = passRealCR  (leptons, jets, m_met, CR); break;
    case sf::CR_SideLow  : passCR = passRealCR  (leptons, jets, m_met, CR); break;
    case sf::CR_SideHigh : passCR = passRealCR  (leptons, jets, m_met, CR); break;
    case sf::CR_HF       : passCR = passHFCR    (leptons, jets, m_met, CR); break;
    case sf::CR_HF_high  : passCR = passHFCR    (leptons, jets, m_met, CR); break;
    case sf::CR_Conv     : passCR = passConvCR  (leptons, jets, m_met    ); break;
    case sf::CR_MCConv   : passCR = passMCReg   (leptons, jets, m_met, CR); break;
    case sf::CR_MCQCD    : passCR = passMCReg   (leptons, jets, m_met, CR); break;
    case sf::CR_MCReal   : passCR = passMCReg   (leptons, jets, m_met, CR); break;
        // Remaining signal and control regions
    default          : passCR = passSignalRegion(leptons,jets,m_met,CR); break;
    } // end switch(CR)
    if( passCR ){
      for(size_t ip=0; ip<m_probes.size(); ++ip) fillRatesHistos(m_probes.at(ip), jets, m_met, cr);
    } // if(passCR)
  } // for(cr)
  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// Plotting Method
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::fillRatesHistos(const Lepton* lep, const JetVector& jets,
                                       const Met* met, size_t regionIndex)
{
    struct FillEffHistos {
        bool n_; double w_;
        LeptonType l_; size_t r_; Chan c_;
        FillEffHistos(bool alsoNum, double weight, LeptonType l, size_t r, Chan c)
            : n_(alsoNum), w_(weight), l_(l), r_(r), c_(c) {}
        void operator () (EffObject* eff_array[kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N],
                          double value) {
            if(EffObject* eff = eff_array[l_][r_][c_])         eff->Fill(n_, w_, value);
            if(c_ != Ch_all)
                if(EffObject* eff = eff_array[l_][r_][Ch_all]) eff->Fill(n_, w_, value);
        }
    };

    bool isElOrMu(lep->isEle() || lep->isMu());
    assert(isElOrMu);
    LeptonType lt(lep->isEle() ? kElectron : kMuon);
    bool pass(isSignalLepton(lep, m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC));
    FillEffHistos fillEff(pass, m_evtWeight, lt, regionIndex, static_cast<Chan>(m_ch));    
    fillEff(h_l_pt        , lep->Pt());
    fillEff(h_l_pt_coarse , lep->Pt());
    fillEff(h_l_eta       , fabs(lep->Eta()));
    fillEff(h_l_eta_coarse, fabs(lep->Eta()));
    fillEff(h_metrel      , m_metRel);
    fillEff(h_met         , met->Et);
    fillEff(h_njets       , jets.size());
    fillEff(h_onebin      , 0.0);
    if( nt.evt()->isMC ){ // If the event is MC, save the flavor
        LeptonSource ls(getLeptonSource(lep));
        bool isQcd(ls==LS_HF||ls==LS_LF);
        fillEff(h_flavor, ls);
        if(isQcd) fillEff(h_flavor, LS_QCD);
    }
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
				 sf::Region CR)
{
  // Used to measure fake rate in MC in a region
  // that is close to our signal region
  // * 40 < metrel < 100
  // * exactly two leptons
  // * probes are the fake leptons of specific type
  if( !nt.evt()->isMC )     return false;
  if( leptons.size() != 2 ) return false;
  m_ch = SusySelection::getChan(leptons);
  m_metRel = getMetRel(met,leptons,jets);
  for(uint il=0; il<leptons.size(); ++il){
    Lepton* l=leptons[il];
    bool isQcdLepton(susy::isHFLepton(l) || susy::isLFLepton(l));
    if(CR==CR_MCConv && susy::isConvLepton(l)) m_probes.push_back( l ); // Conversion
    if(CR==CR_MCQCD && isQcdLepton)            m_probes.push_back( l ); // QCD
    if(CR==CR_MCReal && susy::isRealLepton(l)) m_probes.push_back( l ); // Real
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
                                        sf::Region CR)
{
  if( leptons.size() != 2 ) return false;
  bool passSR = false;
  bool computeWhssPass(CR==susy::fake::CR_SRWHSS || CR==susy::fake::CR_CR8lpt
                       || CR==susy::fake::CR_CR8ee || CR==susy::fake::CR_CR8mm
                       || CR==susy::fake::CR_CR8mmMtww || CR==susy::fake::CR_CR8mmHt);
  SsPassFlags whssFlags(computeWhssPass ? passWhSS(leptons,jets,met) : SsPassFlags());
  susy::wh::Chan ch(SusySelection::getChan(leptons));
  m_ch = SusySelection::getChan(leptons);
  bool isEe(susy::wh::Ch_ee==ch), isMm(susy::wh::Ch_mm==ch);
  switch(CR) {
  case susy::fake::CR_SSInc        : passSR = susy::sameSign(leptons);            break;
  case susy::fake::CR_SRWHSS       : passSR =  whssFlags.metrel;                  break;
  case susy::fake::CR_CR8lpt       : passSR =  whssFlags.lepPt;                   break;
  case susy::fake::CR_CR8ee        : passSR = (whssFlags.lepPt   && isEe);        break;
  case susy::fake::CR_CR8mm        : passSR = (whssFlags.lepPt   && isMm);        break;
  case susy::fake::CR_CR8mmMtww    : passSR = (whssFlags.mtllmet && isMm);        break;
  case susy::fake::CR_CR8mmHt      : passSR = (whssFlags.ht      && isMm);        break;
  case susy::fake::CR_SsEwk        : passSR = SusySelection::passEwkSs     (leptons, jets, met); break;
  case susy::fake::CR_SsEwkLoose   : passSR = SusySelection::passEwkSsLoose(leptons, jets, met); break;
  default: cout<<"invalid ControlRegion "<<CR<<endl;
  }
  for(uint i=0; i<leptons.size(); ++i) m_probes.push_back( leptons[i]);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);
  return passSR;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Real Control Region: Z
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passRealCR(const LeptonVector &leptons,
                                  const JetVector &jets,
                                  const Met* met,
                                  sf::Region CR)
{
  // Real CRA:
  // * Require Exactly two baseline leptons SF
  // * Require single lepton triggers
  // * Require fall in Z mass window
  // * Require at least one tight lepton for tag

  // Same flavor dilepton
  if( leptons.size() != 2 )  return false;
  if( !susy::sameFlavor(leptons) ) return false;
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
                                sf::Region CR)
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
Susy::Lepton* muon_ptr2lepton_ptr(Susy::Muon *m) { return static_cast<Susy::Lepton*>(m); }
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
  LeptonVector tempL(preMuons.size(), NULL);
  std::transform(preMuons.begin(), preMuons.end(), tempL.begin(), muon_ptr2lepton_ptr);
  bool passTrigger(SusySelection::passTrig2LwithMatch(tempL, m_trigObj, m_met->Et, nt.evt()));
  if( !passTrigger )                         return false; // Pt and trigger requirement
  if( met->Et > 50 )                         return false; // Met < 50
  if( (*preMuons[0]+*preMuons[1]).M() < 20 ) return false; // M(mu,mu) > 20
  if( Mt((Lepton*) preElecs[0], met) > 40)   return false; // Mt(elec, met) < 40
  m_probes.push_back( (Lepton*) preElecs[0] );
  m_metRel = getMetRel(met, leptons, jets); // Not useful quantity for this region.
  float eff = nt.evt()->isMC ? m_probes[0]->effSF : 1.0;
  m_evtWeight = getEvtWeight(tempL,true,true) * eff;
  return true;
}
//----------------------------------------------------------
sf::LeptonSource MeasureFakeRate2::getLeptonSource(const Lepton* l)
{
  if( isRealLepton(l) ) return LS_Real;
  if( susy::isHFLepton(l) )   return LS_HF;
  if( susy::isLFLepton(l) )   return LS_LF;
  if( susy::isConvLepton(l) ) return LS_Conv;
  return LS_Unk;
}
//----------------------------------------------------------
const std::vector<susy::fake::Region> MeasureFakeRate2::allRegions() const
{
    vector<sf::Region> allRegions;
    allRegions.insert(allRegions.end(), m_controlRegions.begin(), m_controlRegions.end());
    allRegions.insert(allRegions.end(), m_signalRegions.begin(), m_signalRegions.end());
    return allRegions;
}
//----------------------------------------------------------
const std::string MeasureFakeRate2::LeptonType2str(const LeptonType l)
{
    string lname("unknown");
    switch(l) {
    case kElectron : lname = "elec"; break;
    case kMuon     : lname = "muon"; break;
    default        : ;
    }
    return lname;
}
//----------------------------------------------------------
