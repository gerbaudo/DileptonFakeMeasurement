#include "SusyTest0/MeasureFakeRate2.h"
#include "SusyTest0/criteria.h"
#include "SusyTest0/FakeBinnings.h"
#include "SusyTest0/SusySelection.h" // passEwkSs*
#include "SusyTest0/kinematic.h"
#include "SusyTest0/utils.h"

#include "LeptonTruthTools/RecoTruthMatch.h"

#include <algorithm> // transform
#include <cassert>

using namespace susy::fake;
namespace sf = susy::fake;
using namespace susy::wh;
namespace swk = susy::wh::kin;


const sf::Region controlRegions[] = {
    sf::CR_Real, sf::CR_SideLow, sf::CR_SideHigh, sf::CR_HF, sf::CR_HF_high,
    sf::CR_Conv,
    sf::CR_MCConv, sf::CR_MCQCD,
    sf::CR_MCReal
};
const size_t nControlRegions = sizeof(controlRegions)/sizeof(controlRegions[0]);
const sf::Region signalRegions[] = {
  sf::CR_SSInc,
  sf::CR_SSInc1j,
  sf::CR_SRWHSS,
  sf::CR_CR8lpt,
  sf::CR_CR8ee,
  sf::CR_CR8mm,
  sf::CR_CR8mmMtww,
  sf::CR_CR8mmHt,
  sf::CR_CR9lpt,
  sf::CR_SsEwk,
  sf::CR_SsEwkLoose,
  sf::CR_SsEwkLea,
  sf::CR_WHZVfake1jee,
  sf::CR_WHZVfake2jee,
  sf::CR_WHZVfake1jem,
  sf::CR_WHZVfake2jem,
  sf::CR_WHfake1jem,
  sf::CR_WHfake2jem,
  sf::CR_WHZV1jmm,
  sf::CR_WHZV2jmm,
  sf::CR_WHfake1jmm,
  sf::CR_WHfake2jmm,

  sf::CR_WHZVfake1j,
  sf::CR_WHZVfake2j,
  sf::CR_WHfake1j,
  sf::CR_WHfake2j,
  sf::CR_WHZV1j,
  sf::CR_WHZV2j,

  sf::CR_SRWH1j,
  sf::CR_SRWH2j
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
  m_ch(0),
  m_ET(ET_Unknown)
{
  resetCounters();
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
  SusySelection::Begin(0);
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
        h_l_pt_true    [il][icr][ich] = new TH1F     ((bn+"l_pt_true").c_str(), "", nFakePtbins, FakePtbins);
        h_l_pt_fake    [il][icr][ich] = new TH1F     ((bn+"l_pt_fake").c_str(), "", nFakePtbins, FakePtbins);
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
        TH1F *hlpt = (ls==LS_Real ? h_l_pt_true[lt][regionIndex][m_ch] : h_l_pt_fake[lt][regionIndex][m_ch]);
        if(hlpt) hlpt->Fill(lep->Pt(), m_evtWeight);
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
//  if(m_metRel<40.0) return false;
  if(jets.size()<1) return false;
//  if(jets.size()<1 || m_metRel<20.0) return false;
  for(uint il=0; il<leptons.size(); ++il){
    Lepton* l=leptons[il];
    bool isQcdLepton(susy::isHFLepton(l) || susy::isLFLepton(l));
    uint dsid(nt.evt()->mcChannel);
    if(CR==CR_MCConv && susy::isConvLepton(l)) m_probes.push_back( l ); // Conversion
    if(CR==CR_MCQCD && isQcdLepton)            m_probes.push_back( l ); // QCD
    if(CR==CR_MCReal && isRealLepton(l, dsid)) m_probes.push_back( l ); // Real
  }
  m_evtWeight = MeasureFakeRate2::getEvtWeight(leptons);
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
  bool passSR = false;
  // all the signal regions below require 2L SS
  if( leptons.size() != 2 ) return passSR;
  if(!susy::sameSign(leptons)) return passSR;
  SsPassFlags whssFlags = MeasureFakeRate2::passWhSS(leptons,jets,met); // NB this is not the method from SusySelection. Refactor
  susy::wh::Chan ch(SusySelection::getChan(leptons));
  m_ch = SusySelection::getChan(leptons);
  bool isEe(susy::wh::Ch_ee==ch), isMm(susy::wh::Ch_mm==ch), isOf(!isEe && !isMm);
  bool is1j(jets.size()==1), is2j(jets.size()>1);
  LeptonVector anyLeptons(getAnyElOrMu(nt));
  LeptonVector lowPtLep(subtract_vector(anyLeptons, m_baseLeptons));
  /*const*/ swk::DilepVars v(swk::compute2lVars(leptons, met, jets));
  v.l3veto = whssFlags.veto3rdL = SusySelection::passThirdLeptonVeto(leptons[0], leptons[1], lowPtLep, m_debugThisEvent);
  if(CR==sf::CR_SSInc) passSR = susy::sameSign(leptons);
  else {
      bool passCommonCriteria (whssFlags.sameSign && whssFlags.tauVeto && whssFlags.fjveto && whssFlags.bjveto && whssFlags.ge1j);
      if(passCommonCriteria) {
          switch(CR) {
          case sf::CR_SSInc1j      : passSR =  whssFlags.ge1j;                    break;
          case sf::CR_SRWHSS       : passSR =  whssFlags.metrel;                  break;
          case sf::CR_CR8lpt       : passSR =  whssFlags.lepPt;                   break;
          case sf::CR_CR8ee        : passSR = (whssFlags.lepPt   && isEe);        break;
          case sf::CR_CR8mm        : passSR = (whssFlags.lepPt   && isMm);        break;
          case sf::CR_CR8mmMtww    : passSR = (whssFlags.mtllmet && isMm);        break;
          case sf::CR_CR8mmHt      : passSR = (whssFlags.ht      && isMm);        break;
          case sf::CR_CR9lpt       : passSR = (isOf && whssFlags.lepPt);          break;
          case sf::CR_SsEwk        : passSR = SusySelection::passEwkSs     (leptons, jets, met); break;
          case sf::CR_SsEwkLoose   : passSR = SusySelection::passEwkSsLoose(leptons, jets, met); break;
          case sf::CR_SsEwkLea     : passSR = SusySelection::passEwkSsLea  (leptons, jets, met); break;
          case sf::CR_WHZVfake1jee : passSR = (isEe && is1j && SusySelection::passCrWhZVfakeEe(v)); break;
          case sf::CR_WHZVfake2jee : passSR = (isEe && is2j && SusySelection::passCrWhZVfakeEe(v)); break;
          case sf::CR_WHZVfake1jem : passSR = (isOf && is1j && SusySelection::passCrWhZVfakeEm(v)); break;
          case sf::CR_WHZVfake2jem : passSR = (isOf && is2j && SusySelection::passCrWhZVfakeEm(v)); break;
          case sf::CR_WHfake1jem   : passSR = (isOf && is1j && SusySelection::passCrWhfakeEm(v)); break;
          case sf::CR_WHfake2jem   : passSR = (isOf && is2j && SusySelection::passCrWhfakeEm(v)); break;
          case sf::CR_WHZV1jmm     : passSR = (isMm && is1j && SusySelection::passCrWhZVMm(v)); break;
          case sf::CR_WHZV2jmm     : passSR = (isMm && is2j && SusySelection::passCrWhZVMm(v)); break;
          case sf::CR_WHfake1jmm   : passSR = (isMm && is1j && SusySelection::passCrWhfakeMm(v)); break;
          case sf::CR_WHfake2jmm   : passSR = (isMm && is2j && SusySelection::passCrWhfakeMm(v)); break;
              // trying to unify the ee/em/mm channels
          case sf::CR_WHZVfake1j   : passSR = (is1j && passCrWhZVfake(v)); break;
          case sf::CR_WHZVfake2j   : passSR = (is2j && passCrWhZVfake(v)); break;
          case sf::CR_WHfake1j     : passSR = (is1j && passCrWhfake  (v)); break;
          case sf::CR_WHfake2j     : passSR = (is2j && passCrWhfake  (v)); break;
          case sf::CR_WHZV1j       : passSR = (is1j && passCrWhZV    (v)); break;
          case sf::CR_WHZV2j       : passSR = (is2j && passCrWhZV    (v)); break;

          case sf::CR_SRWH1j       : passSR = (is1j && passSrWh1j    (v)); break;
          case sf::CR_SRWH2j       : passSR = (is2j && passSrWh2j    (v)); break;

          default: cout<<"invalid ControlRegion "<<CR<<endl;
          } // switch
      } // passCommonCriteria
  }
  for(uint i=0; i<leptons.size(); ++i) m_probes.push_back( leptons[i]);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);
  return passSR;
}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::selectEvent(bool count)
{
  if(m_dbg) cout << "MeasureFakeRate2::selectEvent" << endl;
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
  if(!passDeadRegions(m_preJets, m_met, nt.evt()->run, nt.evt()->isMC)) return false;
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
  if( !SusySelection::passHfor(nt) ) return false;
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
SsPassFlags MeasureFakeRate2::passWhSS(const LeptonVector& leptons, const JetVector& jets, const Met* met)
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
  if(susy::sameSign(leptons)) { increment(n_pass_CRWHSS2lss  [m_ET], lsf, bsf); f.sameSign=true;} else return f;
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
float MeasureFakeRate2::getEvtWeight(const LeptonVector& leptons, bool includeBTag, bool includeTrig,
				  bool doMediumpp)
{
  if( !nt.evt()->isMC ) return 1.;
  uint nl = leptons.size();
  float weight = 1;
  weight = getEventWeight(LUMI_A_L, true); // lumi, xs, sumw, pileup
  // Trigger
  float trigW = 1;
  if(!m_useMCTrig && includeTrig){
      trigW  = (nl == 2 ?
                m_trigObj->getTriggerWeight(leptons,
                                            nt.evt()->isMC,
                                            m_met->Et,
                                            m_signalJets2Lep.size(),
                                            nt.evt()->nVtx,
                                            NtSys_NOM)
                : 1.0);
    if(trigW != trigW){ cout<<"\tTrigger weight: "<<trigW<<endl; trigW =0; }// deal with NaN
    if(trigW < 0) trigW = 0;
  }
  float effW   = nl > 0 ? leptons[0]->effSF : 1.;
  effW        *=  nl > 1 ? leptons[1]->effSF : 1.;
  // btag, if included
  float bTag   =  includeBTag ? getBTagWeight(nt.evt()) : 1.;
  return weight * trigW * effW * bTag;
}
//----------------------------------------------------------
float MeasureFakeRate2::getBTagWeight(const Event* evt)
{
  JetVector tempJets;
  for(uint ij=0; ij<m_baseJets.size(); ++ij){
    Jet* jet = m_baseJets.at(ij);
    if( !(jet->Pt() > 20 && fabs(jet->detEta) < JET_ETA_CUT_2L) ) continue;
    tempJets.push_back(jet);
  }
  return bTagSF(evt, tempJets, evt->mcChannel, BTag_NOM);
}
//----------------------------------------------------------
void MeasureFakeRate2::dumpEventCounters()
{
  string v_WT[] = {"Raw","Event","Pileup","Pileup A-B3",
                   "LeptonSF","btagSF","TrigSF","All A-B3", "All A-E"};

  for(int w=0; w<WT_N; ++w){
    cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;
    cout << "SusySelection Event counts for weight: " << v_WT[w] << endl;
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
//----------------------------------------------------------
void MeasureFakeRate2::printCounter(string cut, float counter[ET_N][WT_N], int weight)
{
  cout << cut;
  for(int i=0; i<ET_N-2; ++i)
    cout << "\t" << Form("%10.3f",counter[i][weight]);
  cout << endl;
}
//----------------------------------------------------------
void MeasureFakeRate2::increment(float flag[], bool includeLepSF, bool includeBtag)
{
    flag[WT_Raw]   += 1.0;
    flag[WT_Evt]   += nt.evt()->w;
    flag[WT_PU]    += nt.evt()->w * nt.evt()->wPileup;
    flag[WT_PU1fb] += nt.evt()->w * nt.evt()->wPileupAB3;
    flag[WT_LSF]   += (includeLepSF ? 
                       nt.evt()->w * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF :
                       nt.evt()->w);
    float btag = includeBtag ? getBTagWeight(nt.evt()) : 1.0;
    flag[WT_Btag]  += nt.evt()->w * btag;
    
    float trig = m_baseLeptons.size() == 2 && nt.evt()->isMC && !m_useMCTrig ? 
        m_trigObj->getTriggerWeight(m_baseLeptons,
                                    nt.evt()->isMC,
                                    m_met->Et,
                                    m_signalJets2Lep.size(),
                                    nt.evt()->nVtx,
                                    NtSys_NOM) : 1;
    //cout<<"\tTrigger weight: "<<trig<<endl;
    flag[WT_Trig] += trig * nt.evt()->w;
    bool useSumwMap(true);
    
    float allAE = SusyNtAna::getEventWeight(LUMI_A_L, useSumwMap) * btag * trig;
    allAE = includeLepSF ? allAE * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF : allAE;
    flag[WT_AllAE] += allAE;
}
//----------------------------------------------------------
sf::LeptonSource MeasureFakeRate2::getLeptonSource(const Lepton* l)
{
  uint dsid(nt.evt()->mcChannel);
  if( isRealLepton(l, dsid) ) return LS_Real;
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
bool MeasureFakeRate2::isRealLepton(const Lepton* lep, uint dsid)
{
  // Updated way of handling real and fake leptons using LeptonTruthTools
  // Need to handle new g2ww -- Assume all real for now
  if( dsid == 169471 || dsid == 169472 || dsid == 169473 || dsid == 169474 ||
      dsid == 169475 || dsid == 169476 || dsid == 169477 || dsid == 169478 ||
      dsid == 169479)
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
  if( (dsid>=126892 && dsid<=126895) || (dsid>=147770 && dsid<=147772) ||
      (dsid>=147774 && dsid<=147776)){
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
//----------------------------------------------------------
void MeasureFakeRate2::resetCounters()
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
      //
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
}
