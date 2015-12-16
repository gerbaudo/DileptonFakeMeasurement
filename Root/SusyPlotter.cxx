#include <cassert>

#include "DileptonFakeMeasurement/SusyPlotter.h"

#include "TLorentzVector.h"

#include "SusyNtuple/SusyDefs.h"
#include "DileptonFakeMeasurement/criteria.h"
#include "DileptonFakeMeasurement/kinematic.h"
#include "DileptonFakeMeasurement/utils.h"
#include "SusyMatrixMethod/DiLeptonMatrixMethod.h"

using namespace std;
using namespace Susy;
using namespace susy::wh;
namespace swh = susy::wh;
namespace swk = susy::wh::kin;

// Histogram bins
const float varptbins[] = {0,10,20,30,40,50,70,100,150,200,250};
const int    varnptbins = 10;

const float ptmin    = 0;
const float ptmax    = 250;
const int  nptbins   = 25;

const float etamin   = 0;
const float etamax   = 5;
const int   netabins = 25;

const float dphimin   = 0;
const float dphimax   = 3.2;
const int   ndphibins = 15;

const float drmin   = 0;
const float drmax   = 5;
const int   ndrbins = 25;

const float massmin  = 0;
const float massmax  = 300;
const int  nmassbins = 30;

const float massminfine  = 0;
const float massmaxfine  = 300;
const int  nmassbinsfine = 60;

const float massminfiner  = 0;
const float massmaxfiner  = 300;
const int  nmassbinsfiner = 100;

const float massminj  = 0;
const float massmaxj  = 1000;
const int  nmassbinsj = 100;

const int   njetbins = 7;
const float njetmin = -0.5;
const float njetmax = njetbins - 0.5;

const float Typemin = -0.5;
const float Typemax = 22.5;
const int     nType = 23;

const float Originmin = -0.5;
const float Originmax = 42.5;
const int     nOrigin = 43;
//-----------------------------------------
SusyPlotter::SusyPlotter() :
  SusySelection(),
  m_histFileName("susyPlotterOut.root"),
  m_histFile(0),
  m_doFake(false)
{
}
//-----------------------------------------
void SusyPlotter::Begin(TTree* /*tree*/)
{
  SusySelection::Begin(0);
  if(m_dbg) cout << "SusyPlotter::Begin" << endl;
  setSysts();
  initHistos();
}
//-----------------------------------------
Bool_t SusyPlotter::Process(Long64_t entry)
{
  m_printer.countAndPrint(cout);
  GetEntry(entry);
  clearObjects();
  cacheStaticWeightComponents();
  increment(n_readin, m_weightComponents);
  bool removeLepsFromIso(false);
  const SusyNtSys sys = NtSys_NOM; // for now we don't do any sys loop; see also setSysts
  selectObjects(sys, removeLepsFromIso, TauID_medium);
  if(!selectEvent())    return kTRUE;
  const Met*          m = m_met;
  const JetVector&    j = m_signalJets2Lep;
  const JetVector&   bj = m_baseJets;     // DG don't know why, but we use these for the btag w
  const LeptonVector& l = m_signalLeptons;
  LeptonVector&     ncl = m_signalLeptons; // non-const leptons: can be modified by qflip
  Met ncmet(*m_met); // non-const met
  const TauVector&    t = m_signalTaus;
  if(l.size()>1) computeNonStaticWeightComponents(l, bj); else return false;
  bool allowQflip(true);
  SsPassFlags ssf(SusySelection::passSrSs(WH_SRSS1, ncl, t, j, m, allowQflip));
  if(!ssf.passCommonCriteria()) return false;
  float weight(m_weightComponents.product());
  const DiLepEvtType ll(getDiLepEvtType(l)), ee(ET_ee), mm(ET_mm);
  bool sameFlav(ll==ee||ll==mm);
  if(ssf.passLpt()) {
    swh::Region pr = (sameFlav ? swh::PR_CR8lpt : swh::PR_CR9lpt);
    fillHistos(ncl, j, m, weight, pr, sys);
    if     (ll==ee && ssf.zllVeto) fillHistos(ncl, j, m, weight, PR_CR8ee, sys);
    else if(ll==mm) {
      bool passMinMet(m->Et > 40.0);
      if(passMinMet ) fillHistos(ncl, j, m, weight, swh::PR_CR8mm,     sys);
      if(ssf.mtllmet) fillHistos(ncl, j, m, weight, swh::PR_CR8mmMtww, sys);
      if(ssf.ht     ) fillHistos(ncl, j, m, weight, swh::PR_CR8mmHt,   sys);
    } // end if(mm)
  } // end passLpt
  bool passEwkSs     (SusySelection::passEwkSs     (ncl,j,m));
  bool passEwkSsLoose(SusySelection::passEwkSsLoose(ncl,j,m));
  bool passEwkSsLea  (SusySelection::passEwkSsLea  (ncl,j,m));
  if(passEwkSs)      fillHistos(ncl, j, m, weight, swh::PR_SsEwk,     sys);
  if(passEwkSsLoose) fillHistos(ncl, j, m, weight, swh::PR_SsEwkLoose,sys);
  if(passEwkSsLea)   fillHistos(ncl, j, m, weight, swh::PR_SsEwkLea,  sys);
  bool isEe(ll==ee), isMm(ll==mm), isOf(!isEe && !isMm);
  bool is1j(j.size()==1), is2j(j.size()>1);
  LeptonVector anyLeptons(getAnyElOrMu(nt));
  LeptonVector lowPtLep(subtract_vector(anyLeptons, m_baseLeptons));
  /*const*/ swk::DilepVars v(swk::compute2lVars(ncl, m, j));
  v.l3veto = ssf.veto3rdL; // already computed in passSrSs

  if(ssf.sameSign && ssf.ge1j)                           fillHistos(ncl, j, m, weight, swh::PR_CRSsInc1j, sys);
  if(isEe && is1j && SusySelection::passCrWhZVfakeEe(v)) fillHistos(ncl, j, m, weight, swh::CrZVfake1jee, sys);
  if(isEe && is2j && SusySelection::passCrWhZVfakeEe(v)) fillHistos(ncl, j, m, weight, swh::CrZVfake2jee, sys);
  if(isOf && is1j && SusySelection::passCrWhZVfakeEm(v)) fillHistos(ncl, j, m, weight, swh::CrZVfake1jem, sys);
  if(isOf && is2j && SusySelection::passCrWhZVfakeEm(v)) fillHistos(ncl, j, m, weight, swh::CrZVfake2jem, sys);
  if(isOf && is1j && SusySelection::passCrWhfakeEm  (v)) fillHistos(ncl, j, m, weight, swh::Crfake1jem  , sys);
  if(isOf && is2j && SusySelection::passCrWhfakeEm  (v)) fillHistos(ncl, j, m, weight, swh::Crfake2jem  , sys);
  if(isMm && is1j && SusySelection::passCrWhZVMm    (v)) fillHistos(ncl, j, m, weight, swh::CrZV1jmm    , sys);
  if(isMm && is2j && SusySelection::passCrWhZVMm    (v)) fillHistos(ncl, j, m, weight, swh::CrZV2jmm    , sys);
  if(isMm && is1j && SusySelection::passCrWhfakeMm  (v)) fillHistos(ncl, j, m, weight, swh::Crfake1jmm  , sys);
  if(isMm && is2j && SusySelection::passCrWhfakeMm  (v)) fillHistos(ncl, j, m, weight, swh::Crfake2jmm  , sys);

  if(is1j && SusySelection::passAndIncrementCrWhZVfake(v)) fillHistos(ncl, j, m, weight, swh::CrZVfake1j , sys);
  if(is2j && SusySelection::passAndIncrementCrWhZVfake(v)) fillHistos(ncl, j, m, weight, swh::CrZVfake2j , sys);
  if(is1j && SusySelection::passAndIncrementCrWhfake  (v)) fillHistos(ncl, j, m, weight, swh::Crfake1j   , sys);
  if(is2j && SusySelection::passAndIncrementCrWhfake  (v)) fillHistos(ncl, j, m, weight, swh::Crfake2j   , sys);
  if(is1j && SusySelection::passAndIncrementCrWhZV    (v)) fillHistos(ncl, j, m, weight, swh::CrZV1j     , sys);
  if(is2j && SusySelection::passAndIncrementCrWhZV    (v)) fillHistos(ncl, j, m, weight, swh::CrZV2j     , sys);

  if(is1j && SusySelection::passSrWh1j(v)) fillHistos(ncl, j, m, weight, swh::SrWh1j , sys);
  if(is2j && SusySelection::passSrWh2j(v)) fillHistos(ncl, j, m, weight, swh::SrWh2j , sys);


  Lepton *l0 = ncl[0], *l1 = ncl[1];
  //bool isSuspiciousEvent(isOf && SusySelection::passCrWhfake(v) && fabs(l0->Eta())>1.80);
  bool isEventToPrint = (SusySelection::passCrWhZVfake(v) ||
                         SusySelection::passCrWhfake(v) ||
                         SusySelection::passCrWhZV(v) ||
                         SusySelection::passSrWh1j(v) || SusySelection::passSrWh2j(v));
  if(isEventToPrint) {
      unsigned int run(nt.evt()->run), event(nt.evt()->event);
      cout<<"fake suspicious event run "<<run<<" event "<<event
          <<" "<<(isEe?"ee":isMm?"mm":"em")
          <<" "<<(is1j?"1j":is2j?"2j":"nj")<<endl;
      nt.evt()->print();
      l0->print();
      l1->print();
      m->print();
      for(size_t iJet=0; iJet<j.size(); ++iJet)
          j[iJet]->print();
  }

  return kTRUE;
}
//-----------------------------------------
void SusyPlotter::Terminate()
{
  SusySelection::Terminate();
  if(m_dbg) cout << "SusyPlotter::Terminate" << endl;
  dumpEventCounters();
  // Save the output
  m_histFile->Write();
  m_histFile->Close();
}
//----------------------------------------------------------
SusyPlotter& SusyPlotter::setOutputFilename(const std::string &name)
{
  bool invalidFilename(name.size()<1 || string::npos==name.find(".root"));
  if(invalidFilename)
    cout<<"Warning! SusyPlotter::setOutputFilename('"<<name<<"')"<<" invalid filename."<<endl
        <<"\t using default value '"<<m_histFile<<"'"<<endl;
  else m_histFileName = name;
  return *this;
}
//-----------------------------------------
void SusyPlotter::fillHistos(const LeptonVector& leps, const JetVector &jets,
                             const Met* met, const float weight,
                             size_t regionIndex, uint sys)
{
  if(m_dbg) cout << "SusyPlotter::fillHistos" << endl;
  if( leps.size() != 2 ) return;
  susy::wh::Chan ch = SusySelection::getChan(leps);
  const Lepton* l0 = leps[0];
  const Lepton* l1 = leps[1];
  const size_t &ri = regionIndex;
  assert(l0);
  assert(l1);
  #define FILL(h, var)					\
    do{								\
      float max   = h[ch][ri][sys]->GetXaxis()->GetXmax();	\
      float xfill = var > max ? max - 1e-4 : var;		\
      h[ch][ri][sys]->Fill(xfill,weight);			\
      h[Ch_all][ri][sys]->Fill(xfill,weight);			\
    }while(0)

  #define FILL2(h, varx, vary)						\
    do{									\
      float maxx   = h[ch][ri][sys]->GetXaxis()->GetXmax();	\
      float xfill = varx > maxx ? maxx - 1e-4 : varx;		\
      float maxy   = h[ch][ri][sys]->GetYaxis()->GetXmax();	\
      float yfill = vary > maxy ? maxy - 1e-4 : vary;		\
      h[ch][ri][sys]->Fill(xfill,yfill,weight);			\
      h[Ch_all][ri][sys]->Fill(xfill,yfill,weight);			\
    }while(0)

  const TLorentzVector mlv = met->lv();
  const TLorentzVector ll  = *l0 + *l1;

  FILL(h_onebin, 0.);
  FILL(h_l0_pt, l0->Pt());
  FILL(h_l1_pt, l1->Pt());
  FILL(h_ll_pt, ll.Pt());
  if(l0->isEle()) FILL(h_e_pt, l0->Pt());
  else            FILL(h_m_pt, l0->Pt());
  if(l1->isEle()) FILL(h_e_pt, l1->Pt());
  else            FILL(h_m_pt, l1->Pt());

  FILL(h_l0_eta, l0->Eta());
  FILL(h_l1_eta, l1->Eta());
  if(l0->isEle()) FILL(h_e_eta, fabs(l0->Eta()));
  else            FILL(h_m_eta, fabs(l0->Eta()));
  if(l1->isEle()) FILL(h_e_eta, fabs(l1->Eta()));
  else            FILL(h_m_eta, fabs(l1->Eta()));

  FILL(h_ll_M, ll.M());

  float metrel = getMetRel(met, leps, jets);
  FILL(h_met, met->Et);
  FILL(h_metrel, metrel);

  int nJ = jets.size() ; //> 4 ? 4 : jets.size();
  int nBaseJets = m_baseJets.size();
  int nbJ = numberOfCBJets(jets);
  int nfJ = numberOfFJets(jets);

  FILL(h_njets,     nJ);
  FILL(h_nbasejets, nBaseJets);
  FILL(h_nbjets,    nbJ);
  FILL(h_nfjets,    nfJ);
  if(l0->q > 0) FILL(h_njets_pos,nJ);
  else          FILL(h_njets_neg, nJ);

  FILL(h_l_type, l0->mcType);
  FILL(h_l_origin, l0->mcOrigin);
  FILL(h_l_type, l1->mcType);
  FILL(h_l_origin, l1->mcOrigin);

  FILL(h_sumQ, l0->q + l1->q);
  float mt_met_ll = susy::transverseMass(ll, mlv);
  float mt2 = susy::computeMt2(*l0, *l1, mlv);

  FILL(h_met_ll_M, (*l0 + *l1 + mlv).M());
  FILL(h_met_ll_Mt, mt_met_ll);
  FILL(h_mtautau_l0l1met, susy::mZTauTau(*l0, *l1, mlv));
  FILL(h_mt_ll_met, susy::transverseMass(ll, mlv));
  FILL(h_met_ll_Mt2, mt2);

  FILL(h_dR_l0_l1, fabs(l0->DeltaR(*l1)));
  FILL(h_dPhi_l0_l1, fabs(l0->DeltaPhi(*l1)));
  FILL(h_dPhi_met_l0, fabs(met->lv().DeltaPhi(*l0)));
  FILL(h_dPhi_met_l1, fabs(met->lv().DeltaPhi(*l1)));
  FILL(h_dPhi_met_ll, fabs(met->lv().DeltaPhi(*l0+*l1)));
  FILL(h_l0_qeta, (l0->q > 0. ? +1.0 : -1.0) * l0->Eta());
  FILL(h_l1_qeta, (l1->q > 0. ? +1.0 : -1.0) * l1->Eta());

  float mt_l0 = Mt(*l0, mlv);
  float mt_l1 = Mt(*l1, mlv);
  FILL(h_mt_l0_met, mt_l0);
  FILL(h_mt_l1_met, mt_l1);
  FILL(h_mt_l_met_min, (mt_l0 < mt_l1 ? mt_l0 : mt_l1));
  bool topTag = passTopTag(leps, jets, met);
  FILL(h_mct_top_tag, float(topTag));

  if(nJ>=2){
    const Jet &j0 = *jets.at(0);
    const Jet &j1 = *jets.at(1);
    const TLorentzVector jj = j0 + j1;
    FILL(h_j0_pt, j0.Pt());
    FILL(h_j1_pt, j1.Pt());
    FILL(h_j0_eta, j0.Eta());
    FILL(h_j1_eta, j1.Eta());
    FILL(h_jj_deta, fabs(j0.Eta()-j1.Eta()));
    FILL(h_jj_drap, fabs(j0.Rapidity()-j1.Rapidity()));
    FILL(h_jj_M, jj.M());
    FILL(h_dR_ll_jj, fabs(ll.DeltaR(jj)));
    FILL(h_dPhi_ll_jj, fabs(ll.DeltaPhi(jj)));
    FILL(h_dPhi_l0_jj, fabs(l0->DeltaPhi(jj)));
    FILL(h_dPhi_l1_jj, fabs(l1->DeltaPhi(jj)));
    FILL(h_tot_pt, (ll+jj+mlv).Pt());
    FILL(h_sumJ0J1_mv1tag, j0.mv1 + j1.mv1);
    bool oppositeCharge((l0->q * l1->q) < 0.);
    if(oppositeCharge) {
      const TLorentzVector &lPos = (l0->q > 0. ? *l0 : *l1);
      const TLorentzVector &lNeg = (l1->q < 0. ? *l1 : *l0);
      int numNeutrinoSol = 0;
      numNeutrinoSol += susy::numberOfNeutrinoSolutions(lPos, lNeg, j0, j1, mlv);
      numNeutrinoSol += susy::numberOfNeutrinoSolutions(lPos, lNeg, j1, j0, mlv);
      // DG Mar13: should also we consider the combinations with the third jet?
      FILL(h_numNeutrinoSol, numNeutrinoSol);
    } // end if(oppositeCharge)
  } // end if(nJ>=2)
  #undef FILL
  #undef FILL2
}
//-----------------------------------------
void SusyPlotter::setSysts()
{
  if(!m_doFake) {
    m_systs.push_back(NtSys_NOM);  m_systNames.push_back(SusyNtSystNames[NtSys_NOM]);
  } else {
    // DG 2013-10-18: for now we are not doing any syst loop. However,
    // we need to book the histograms for the systematic variations
    // used to determine the uncertainty on the fake estimate (-> m_doFake toggle).
    // This implementation is really confusing and not safe because
    // here we are potentially mixing two enums
    // (SusyMatrixMethod::SYSTEMATIC and SusyNtSys from
    // SusyNtuple/SusyDefs.h). This will be fixed when I move to the
    // 'single-enum' implementation.
    namespace smm = SusyMatrixMethod;
    const std::string *sns = smm::systematic_names;
    m_systs.push_back(smm::SYS_NOM);        m_systNames.push_back(sns[smm::SYS_NOM]);
    m_systs.push_back(smm::SYS_EL_RE_UP);    m_systNames.push_back(sns[smm::SYS_EL_RE_UP]);
    m_systs.push_back(smm::SYS_EL_RE_DOWN);  m_systNames.push_back(sns[smm::SYS_EL_RE_DOWN]);
    m_systs.push_back(smm::SYS_MU_RE_UP);    m_systNames.push_back(sns[smm::SYS_MU_RE_UP]);
    m_systs.push_back(smm::SYS_MU_RE_DOWN);  m_systNames.push_back(sns[smm::SYS_MU_RE_DOWN]);
    m_systs.push_back(smm::SYS_EL_FR_UP);    m_systNames.push_back(sns[smm::SYS_EL_FR_UP]);
    m_systs.push_back(smm::SYS_EL_FR_DOWN);  m_systNames.push_back(sns[smm::SYS_EL_FR_DOWN]);
    m_systs.push_back(smm::SYS_MU_FR_UP);    m_systNames.push_back(sns[smm::SYS_MU_FR_UP]);
    m_systs.push_back(smm::SYS_MU_FR_DOWN);  m_systNames.push_back(sns[smm::SYS_MU_FR_DOWN]);
  }
}
//-----------------------------------------
void SusyPlotter::initHistos()
{
  m_histFile = new TFile(m_histFileName.c_str(), "recreate");
  TH1::SetDefaultSumw2(true);
  m_histFile->cd();
  //m_histFile->mkdir( sysNames[sys].c_str() ) -> cd();
  for(size_t iPR=0; iPR<swh::kNumberOfPlotRegions; ++iPR){   // for(Plot Region)
      string PR(swh::region2str(swh::PlotRegions[iPR]));
    for(uint iCh=0; iCh<susy::wh::Ch_N; ++iCh){ // for(lepton channel)
      string chan = chanNames[iCh];
      for(uint iSys=0; iSys<m_systs.size(); ++iSys){
        string sys = m_systNames.at(iSys);

	// Preprocessor convenience
	// make a histogram by name (leave off the "h_") and binning
        #define NEWHIST(name, xLbl, nbin, min, max)				\
	  do{								\
	    h_ ## name[iCh][iPR][iSys] = new TH1F((PR+"_"+chan+"_"+#name+"_"+sys).c_str(), #name ";" xLbl, nbin, min, max); \
	    h_ ## name[iCh][iPR][iSys]->Sumw2();				\
	  }while(0)
        #define NEWHIST2(name, xLbl, nbin, min, max)				\
	  do{								\
	    h_ ## name[iCh][iPR][iSys] = new TH2F((PR+"_"+chan+"_"+#name+"_"+sys).c_str(), #name ";" xLbl, nbin, min, max, nbin, min, max); \
	    h_ ## name[iCh][iPR][iSys]->Sumw2();				\
	  }while(0)

        #define NEWVARHIST(name, xLbl, nbin, bins)				\
	  do{								\
	    h_ ## name[iCh][iPR][iSys] = new TH1F((PR+"_"+chan+"_"+#name+"_"+sys).c_str(), #name ";" xLbl, nbin, bins); \
	    h_ ## name[iCh][iPR][iSys]->Sumw2();				\
	  }while(0)

         // Pt Plots
	NEWHIST(l0_pt, "l_{0} P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(j0_pt, "j_{0} P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(l1_pt, "l_{1} P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(j1_pt, "j_{1} P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(e_pt, "Electron P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(m_pt, "Muon P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(ll_pt, "ll P_{T}", nptbins, ptmin, ptmax);
	NEWHIST(tot_pt, "ll+jj+met P_{T}", nptbins, ptmin, ptmax);


	// Eta Plots
	NEWHIST(l0_eta, "l_{0} #eta", netabins, etamin, etamax);
	NEWHIST(j0_eta, "j_{0} #eta", netabins, etamin, etamax);
	NEWHIST(l1_eta, "l_{1} #eta", netabins, etamin, etamax);
	NEWHIST(j1_eta, "j_{1} #eta", netabins, etamin, etamax);
	NEWHIST(jj_deta,"#Delta #eta (j,j)", netabins, etamin, etamax);
	NEWHIST(jj_drap,"#Delta y (j,j)", netabins, etamin, etamax);

	NEWHIST(e_eta, "Electron #eta", netabins, etamin, etamax);
	NEWHIST(m_eta, "Muon #eta", netabins, etamin, etamax);

	// Mass plots
	NEWHIST(ll_M,        "m(ll)",           nmassbins, massmin, massmax);
	NEWHIST(llj_M,       "m(llj)",          nmassbins, massmin, massmax);
	NEWHIST(met_j_M,     "m(met,j)",        nmassbins, massmin, massmax);
	NEWHIST(met_ll_M,    "m(met,ll)",       nmassbins, massmin, massmax);
	NEWHIST(met_j_ll_M,  "m(met,j,ll)",     nmassbins, massmin, massmax);
	NEWHIST(met_j_Mt,    "m_{T}(met,j)",    nmassbins, massmin, massmax);
	NEWHIST(met_ll_Mt,   "m_{T}(met,ll)",   nmassbins, massmin, massmax);
	NEWHIST(met_j_ll_Mt, "m_{T}(met,j,ll)", nmassbins, massmin, massmax);
	NEWHIST(met_ll_Mt2,  "m_{T2}(met,l,l)", nmassbins, massmin, massmax);
	NEWHIST(jj_M,        "m(jj)",           nmassbins, massmin, massmax);

	NEWHIST(mt_ll_met,       "m_{T}(ll,met)", nmassbins, massmin, massmax);
	NEWHIST(mt_l0_met,       "m_{T}(l_{0}, met)", nmassbins, massmin, massmax);
	NEWHIST(mt_l1_met,       "m_{T}(l_{1}, met)", nmassbins, massmin, massmax);
	NEWHIST(mt_l_met_min,    "m_{T}^{min}(l, met)", nmassbins, massmin, massmax);
	NEWHIST(mtautau_l0l1met, "m_{#tau#tau}(l0,l1,met)", nmassbins+1, massmin-((massmax-massmin)/nmassbins), massmax); // need a negative bin
	NEWHIST(mct_top_tag,     "m_{CT} top tag", 2, -0.5, +1.5);
	NEWHIST(sumJ0J1_mv1tag,  "MV1(j0) + MV1(j1)", 50, +0.0, +0.5);
	NEWHIST(numNeutrinoSol,  "Number of neutrino solutions (2l+2j)", 9, -0.5, +8.5);


	// Met plots
	NEWHIST(met, "#slash{E}_{T}", nptbins, ptmin, ptmax);
	NEWHIST(metrel, "#slash{E}^{rel}_{T}", nptbins, ptmin, ptmax);

	// njet plots
	NEWHIST(njets, "# jets", njetbins, njetmin, njetmax);
	NEWHIST(nbasejets, "# base jets", njetbins, njetmin, njetmax);
	NEWHIST(njets_pos, "# jets leading pos", njetbins, njetmin, njetmax);
	NEWHIST(njets_neg, "# jets leading neg", njetbins, njetmin, njetmax);
	NEWHIST(nbjets, "# b jets", njetbins, njetmin, njetmax);
	NEWHIST(nfjets, "# f jets", njetbins, njetmin, njetmax);

	// Type and origin

	NEWHIST(l_type, "l_type", nType, Typemin, Typemax);
	NEWHIST(l_origin, "l_origin", nOrigin, Originmin, Originmax);

	// One bin for counting
	NEWHIST(onebin, "onebin", 1, -0.5, 0.5);

	// Sum charge
	NEWHIST(sumQ, "SumQ", 5, -2.5, 2.5);

	NEWHIST(dPhi_llmet_j, "dPhi(llmet,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dR_llmet_j,   "dR(llmet,j)", ndrbins,drmin, drmax);

	NEWHIST(dPhi_met_l0, "dPhi(met,l0)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_met_l1, "dPhi(met,l1)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_met_ll, "dPhi(met,ll)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_met_j,  "dPhi(met,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_ll_j,   "dPhi(ll,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_l0_j,   "dPhi(l0,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_l1_j,   "dPhi(l1,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_l0_l1,  "dPhi(l0,l1)", ndphibins, dphimin, dphimax);

	NEWHIST(dPhi_woSig_llmet_j, "dPhi(llmet,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_met_l0,  "dPhi(met,l0)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_met_l1,  "dPhi(met,l1)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_met_ll,  "dPhi(met,ll)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_met_j,   "dPhi(met,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_ll_j,    "dPhi(ll,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_l0_j,    "dPhi(l0,j)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_woSig_l1_j,    "dPhi(l1,j)", ndphibins, dphimin, dphimax);

	NEWHIST(dPhi_woSig_l0_l1, "dPhi(l0,l1)", ndphibins, dphimin, dphimax);

	NEWHIST(dPhi_ll_jj, "dPhi(ll,jj)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_l0_jj, "dPhi(l0,jj)", ndphibins, dphimin, dphimax);
	NEWHIST(dPhi_l1_jj, "dPhi(l1,jj)", ndphibins, dphimin, dphimax);

	NEWHIST(dR_l0_l1, "dR(l,l)", ndrbins,drmin, drmax);
	NEWHIST(dR_ll_jj, "dR(ll,jj)", ndrbins,drmin, drmax);

	NEWHIST(l0_qeta, "l_{0} sign(q) #eta", 2*netabins, -etamax, etamax);
	NEWHIST(l1_qeta, "l_{1} sign(q) #eta", 2*netabins, -etamax, etamax);


	NEWHIST2(l0_l1_pt, "l0 vs l1 pt", nptbins, ptmin, ptmax);

        #undef NEWHIST
        #undef NEWHIST2
        #undef NEWVARHIST

      }// end loop over systematics
    }// end loop over channels
  }// end loop over Plot regions
}
//-----------------------------------------
