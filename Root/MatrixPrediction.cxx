#include "SusyTest0/MatrixPrediction.h"

#include "SusyNtuple/SusyDefs.h"
#include "SusyTest0/utils.h"
#include "SusyTest0/criteria.h"
#include "SusyTest0/kinematic.h"

#include <iomanip>
#include <sstream>      // std::ostringstream

using namespace std;
using namespace Susy;
using namespace susy::wh; // cleanup
namespace swh = susy::wh;

const float ptmin    = 0;
const float ptmax    = 250;
const int  nptbins   = 25;

const float etamin   = -3;
const float etamax   = 3;
const int   netabins = 30;

const float massmin  = 0;
const float massmax  = 300;
const int  nmassbins = 30;

//----------------------------------------------------------
MatrixPrediction::MatrixPrediction() :
  SusyPlotter(),
  m_matrix(0),
  m_use2dparametrization(false),
  m_allconfigured(false)
{
}
//----------------------------------------------------------
void MatrixPrediction::Begin(TTree* /*tree*/)
{
    if(m_dbg) cout<<"MatrixPrediction::Begin"<<endl;
    SusySelection::Begin(0);
    if(m_writeTuple) {
        // shouldn't we SusyPlotter::toggleNominal() here?
        m_allconfigured = initMatrixTool();
    } else {
        toggleFakeSystematics();
        SusyPlotter::initHistos();
        m_allconfigured = (initMatrixTool() && bookFakeHisto());
    }
}
//----------------------------------------------------------
Bool_t MatrixPrediction::Process(Long64_t entry)
{
  namespace smm = SusyMatrixMethod;
  namespace sf = susy::fake;
  namespace swk = susy::wh::kin;
  namespace sw = susy::wh;
  if(!m_allconfigured) return false;
  m_printer.countAndPrint(cout);
  GetEntry(entry);
  clearObjects();
  cacheStaticWeightComponents();
  increment(n_readin, m_weightComponents);
  bool removeLepsFromIso(false);
  selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium);
  swh::EventFlags eventFlags = computeEventFlags();
  incrementCounters(eventFlags, m_weightComponents);
  if(eventFlags.failAny()) return kTRUE;
  const Met*          m = m_met;
  const JetVector&    j = m_signalJets2Lep;
  const LeptonVector& l = m_baseLeptons;
  LeptonVector&     ncl = m_baseLeptons;
  const TauVector&    t = m_signalTaus;
  if(l.size()<2) return false; // otherwise cannot compute some of the variables
  float metRel = getMetRel(m, l, j);
  DiLepEvtType ll(getDiLepEvtType(l)), ee(ET_ee), mm(ET_mm);
  bool allowQflip(false), passMinMet(m->Et > 40.0);
  bool isEe(ll==ee), isMm(ll==mm), isSf(isEe||isMm), isOf(!isEe && !isMm);
  const VarFlag_t varsFlags(SusySelection::computeSsFlags(ncl, t, j, m, allowQflip));
  const swk::DilepVars &v = varsFlags.first;
  const SsPassFlags &ssf = varsFlags.second;
  m_weightComponents.fake = getFakeWeight(l,sf::CR_SSInc1j, metRel, smm::SYS_NONE); // just for the counters, use generic CR_SRWHSS
  incrementSsCounters(ssf, m_weightComponents);
  if(!ssf.passCommonCriteria()) return false;
  if(m_writeTuple && ssf.lepPt) {
      double weight(getFakeWeight(l, sf::CR_CR8lpt, metRel, smm::SYS_NONE));
      unsigned int run(nt.evt()->run), event(nt.evt()->event);
      LeptonVector anyLep(getAnyElOrMu(nt));
      LeptonVector lowPtLep(subtract_vector(anyLep, m_baseLeptons));
      const Lepton *l0(l[0]), *l1(l[1]);
      const JetVector clJets(SusySelection::filterClJets(m_signalJets2Lep));
      m_tupleMaker.fill(weight, run, event, *l0, *l1, *m, lowPtLep, clJets);
  } else {
      for(uint s = 0; s<m_systs.size(); ++s){
          smm::SYSTEMATIC sys = static_cast<smm::SYSTEMATIC>(m_systs.at(s));
          bool is1j(j.size()==1), is2j(j.size()>1);

          const sf::Region fakeR = sf::CR_SSInc1j;
          if(ssf.sameSign && ssf.ge1j) fillHistos(ncl, j, m, getFakeWeight(l,fakeR,metRel,sys), PR_CRSsInc1j, sys);
          if(isSf && ssf.lepPt    ) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR8lpt,    sys);
          if(isEe && ssf.zllVeto  ) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR8ee,     sys);
          if(isMm && ssf.lepPt
             && passMinMet   )      fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR8mm,     sys);
          if(isMm && ssf.mtllmet  ) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR8mmMtww, sys);
          if(isMm && ssf.ht       ) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR8mmHt,   sys);
          if(isSf && ssf.lepPt    ) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR8lpt,    sys);
          if(isOf && ssf.lepPt    ) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_CR9lpt,    sys);

          if(is1j && SusySelection::passCrWhZVfake(v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::CrZVfake1j , sys);
          if(is2j && SusySelection::passCrWhZVfake(v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::CrZVfake2j , sys);
          if(is1j && SusySelection::passCrWhfake  (v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::Crfake1j   , sys);
          if(is2j && SusySelection::passCrWhfake  (v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::Crfake2j   , sys);
          if(is1j && SusySelection::passCrWhZV    (v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::CrZV1j     , sys);
          if(is2j && SusySelection::passCrWhZV    (v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::CrZV2j     , sys);

          if(is1j && SusySelection::passSrWh1j    (v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::SrWh1j     , sys);
          if(is2j && SusySelection::passSrWh2j    (v)) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), swh::SrWh2j     , sys);

          bool passEwkSs     (SusySelection::passEwkSs     (ncl,j,m));
          bool passEwkSsLoose(SusySelection::passEwkSsLoose(ncl,j,m));
          bool passEwkSsLea  (SusySelection::passEwkSsLea  (ncl,j,m));
          if(passEwkSs)      fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_SsEwk,     sys);
          if(passEwkSsLoose) fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_SsEwkLoose,sys);
          if(passEwkSsLea)   fillHistos(ncl, j, m, getFakeWeight(l,fakeR, metRel,sys), PR_SsEwkLea,  sys);
          if(m_fillHft) {
              bool passWh1j(SusySelection::passSrWh1j(v)), passWh2j(SusySelection::passSrWh2j(v));
              bool storeEvent(passWh1j||passWh2j);
              if(!isHftFillerInitialized()) initHftFiller();
              if(storeEvent) {
                  swh::kin::DilepVars vv(v); // non-const
                  vv.weight = getFakeWeight(l,fakeR, metRel,sys);
                  fillHft(sys, vv);
              }
          }
      } // end for(s)
  }
  return true;
}
//----------------------------------------------------------
void MatrixPrediction::Terminate()
{
    if(m_writeTuple) SusySelection::Terminate();
    else             SusyPlotter::Terminate();
  if(m_dbg) cout << "MatrixPrediction::Terminate" << endl;
  delete m_matrix;
}
//----------------------------------------------------------
bool MatrixPrediction::bookFakeHisto()
{
  if(!m_histFile) {
    cout<<"MatrixPrediction::bookFakeHisto() invalid hist file"<<endl;
    return false;
  }
  m_histFile->cd();  // Histogram file from SusyPlotter
  for(uint iPR=0; iPR<kNumberOfPlotRegions; ++iPR){ // Plot Region
      string PR = swh::region2str(PlotRegions[iPR]);
    for(uint iCh=0; iCh<Ch_N; ++iCh){ // lepton channel
      string chan = chanNames[iCh];
      for(uint iMP=0; iMP<MP_N; ++iMP){  // matrix pairs
        string MP = MPNames[iMP];
        for(uint iWT=0; iWT<WTog_N; ++iWT){ // weight toggle
          string WT = WTNames[iWT];
          for(uint iSYS=0; iSYS<m_systs.size(); ++iSYS){
            string SYS = m_systNames.at(iSYS);
            string base = PR+"_"+chan+"_"+MP+"_"+WT+"_"+SYS;
// Preprocessor convenience: make a histogram by name (leave off the "h_") and binning
#define NEWHIST(name, xLbl, nbin, min, max)                               \
do{                                                                       \
  hf_ ## name[iCh][iPR][iMP][iWT][iSYS]                                   \
    = new TH1F((base+"_"+#name).c_str(), #name ";" xLbl, nbin, min, max); \
  hf_ ## name[iCh][iPR][iMP][iWT][iSYS]->Sumw2();                         \
 }while(0)

            NEWHIST(l0_pt,      "Lepton P_{T}",        nptbins,   ptmin,   ptmax); // Lepton Kin
            NEWHIST(l1_pt,      "Lepton P_{T}",        nptbins,   ptmin,   ptmax);
            NEWHIST(ll_M,       "m(ll)",               nmassbins, massmin, massmax); // Mass
            NEWHIST(met,        "#slash{E}_{T}",       nptbins,   ptmin,   ptmax);  // Met
            NEWHIST(metrel,     "#slash{E}^{rel}_{T}", nptbins,   ptmin,   ptmax);
            NEWHIST(onebin,     "onebin",              1,         -0.5,    0.5);
            NEWHIST(evt_weight, "evt_weight",          5000,      -5,      5);
            NEWHIST(met_l0_Mt,  "met_l0_Mt",           nmassbins, massmin, massmax); // Mt plots
            NEWHIST(met_l1_Mt,  "met_l1_Mt",           nmassbins, massmin, massmax);
            NEWHIST(njets,      "njets",               5,         -0.5,    4.5); // Jet plots
            NEWHIST(bjet_pt,    "bjet_pt",             nptbins,   ptmin,   ptmax);
            NEWHIST(ljet_pt,    "ljet_pt",             nptbins,   ptmin,   ptmax);
#undef NEWHIST
          }// end for(iSYS)
        }// end for(iWT)
      }// end for(iMP)
    }// end for(iCh)
  }// end for(iPR)
  return true;
}
//----------------------------------------------------------
void MatrixPrediction::fillFakeHistos(const LeptonVector &baseLeps, const JetVector &jets,
                                      const Met* met,float weight, size_t regionIndex, uint sys)
{
  const size_t &ri = regionIndex;
  if(m_dbg) cout << "MatrixPrediction::plotFakeHisto" << endl;
  if( baseLeps.size() != 2 ) return;
  int ch = SusySelection::getChan(baseLeps);
  int mp = getMatrixPair(baseLeps);
  const Lepton* l0 = baseLeps[0];
  const Lepton* l1 = baseLeps[1];
#define FILL(h, var)            \
do{                                                                   \
  float max   = h[ch][ri][mp][WT_ON][sys]->GetXaxis()->GetXmax();       \
  float xfill = var > max ? max - 1e-4 : var;                           \
  h[ch][ri][mp][WT_ON][sys]->Fill(xfill,weight);                        \
  h[Ch_all][ri][mp][WT_ON][sys]->Fill(xfill,weight);                    \
  h[ch][ri][mp][WT_OFF][sys]->Fill(xfill,1.0);                          \
  h[Ch_all][ri][mp][WT_OFF][sys]->Fill(xfill,1.0);                      \
  h[ch][ri][MP_ALL][WT_ON][sys]->Fill(xfill,weight);                    \
  h[Ch_all][ri][MP_ALL][WT_ON][sys]->Fill(xfill,weight);                \
  h[ch][ri][MP_ALL][WT_OFF][sys]->Fill(xfill,1.0);                      \
  h[Ch_all][ri][MP_ALL][WT_OFF][sys]->Fill(xfill,1.0);                  \
 }while(0)

  float metrel = getMetRel(met, baseLeps, jets);
  FILL( hf_l0_pt, l0->Pt() );
  FILL( hf_l1_pt, l1->Pt() );
  FILL( hf_ll_M, (*l0 + *l1).M() );
  FILL( hf_met, met->Et );
  FILL( hf_metrel, metrel);
  FILL( hf_onebin, 0. );
  FILL( hf_evt_weight, weight);
  FILL( hf_njets, jets.size() );
  for(uint ij=0; ij<jets.size(); ++ij){
    Jet* jet = jets.at(ij);
    if( isCentralBJet(jet) ) FILL(hf_bjet_pt, jet->Pt());
    if( isCentralLightJet(jet) ) FILL(hf_ljet_pt, jet->Pt());
  }
  #undef FILL
}
//----------------------------------------------------------
float MatrixPrediction::getFakeWeight(const LeptonVector &baseLeps,
                                      susy::fake::Region region,
                                      float metRel,
                                      SusyMatrixMethod::SYSTEMATIC sys)
{
  if(baseLeps.size() != 2) return 0.0;
  uint nVtx = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  float gev2mev(1000.);
  //m_matrix->setDileptonType(baseLeps[0]->isEle(), baseLeps[1]->isEle());
  const Susy::Lepton *l0=baseLeps[0], *l1=baseLeps[1];
  bool l0IsSig(SusyNtTools::isSignalLepton(l0, m_baseElectrons, m_baseMuons, nVtx, isMC));
  bool l1IsSig(SusyNtTools::isSignalLepton(l1, m_baseElectrons, m_baseMuons, nVtx, isMC));
  return m_matrix->getTotalFake(l0IsSig, l0->isEle(), l0->Pt()*gev2mev, l0->Eta(),
                                l1IsSig, l1->isEle(), l1->Pt()*gev2mev, l1->Eta(),
                                region, metRel*gev2mev, sys);
}
//----------------------------------------------------------
float MatrixPrediction::getRFWeight(const LeptonVector &baseLeps,
                                    susy::fake::Region region,
                                    float metRel,
                                    SusyMatrixMethod::SYSTEMATIC sys)
{

  if(baseLeps.size() != 2) return 0.0;
  uint nVtx = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  return m_matrix->getRF( isSignalLepton(baseLeps[0],m_baseElectrons, m_baseMuons,nVtx,isMC),
			  baseLeps[0]->isEle(),
			  baseLeps[0]->Pt() * 1000.,
			  baseLeps[0]->Eta(),
			  isSignalLepton(baseLeps[1],m_baseElectrons, m_baseMuons,nVtx,isMC),
			  baseLeps[1]->isEle(),
			  baseLeps[1]->Pt() * 1000.,
			  baseLeps[1]->Eta(),
			  region,
			  metRel * 1000.,
			  sys);
}
//----------------------------------------------------------
MatrixPair MatrixPrediction::getMatrixPair(const LeptonVector &baseLeps)
{
  if(baseLeps.size() < 2) return MP_N;
  uint nVtx = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  bool l0_tight = isSignalLepton(baseLeps[0],m_baseElectrons, m_baseMuons,nVtx,isMC);
  bool l1_tight = isSignalLepton(baseLeps[1],m_baseElectrons, m_baseMuons,nVtx,isMC);
  if(l0_tight && l1_tight)   return MP_TT;
  if(l0_tight && !l1_tight)  return MP_TL;
  if(!l0_tight && l1_tight)  return MP_LT;
  if(!l0_tight && !l1_tight) return MP_LL;
  return MP_N;
}
//----------------------------------------------------------
MatrixPrediction& MatrixPrediction::setMatrixFilename(const std::string filename)
{
  if(!fileExists(filename))
    cout<<"MatrixPrediction::setMatrixFilename: invalid file '"<<filename<<"'"<<endl
        <<"\t"<<"something will go wrong"<<endl;
  m_matrixFilename = filename;
  return *this;
}
//----------------------------------------------------------
bool MatrixPrediction::initMatrixTool()
{
  // Load the matrix method package
  m_matrix = new SusyMatrixMethod::DiLeptonMatrixMethod();
  SusyMatrixMethod::RATE_PARAM pm = (m_use2dparametrization ? SusyMatrixMethod::PT_ETA : SusyMatrixMethod::PT);
  return m_matrix->configure(m_matrixFilename, pm, pm, pm, pm);
}
//----------------------------------------------------------
std::string MatrixPrediction::dilepDetails(const Susy::Event &event,
                                           const DiLepEvtType &ll,
                                           const LeptonVector &ls)
{
  bool ee(ll==ET_ee), mm(ll==ET_mm);
  const Lepton *l0(ls.size()>0 ? ls[0] : NULL), *l1(ls.size()>1 ? ls[1] : NULL);
  float l0pt(l0 ? l0->Pt() : 0.0), l0eta(l0 ? l1->Eta() : 0.0);
  float l1pt(l1 ? l1->Pt() : 0.0), l1eta(l1 ? l1->Eta() : 0.0);
  std::ostringstream oss;
  oss<<"run "<<event.run
     <<" evt "<<event.event
     <<" "<<(ee?"ee":(mm?"mm":"em"))
     <<" l0: pt="<<l0pt<<" eta="<<l0eta
     <<" l1: pt="<<l1pt<<" eta="<<l1eta;
  return oss.str();
}
//----------------------------------------------------------
std::string MatrixPrediction::eventDetails(bool passSrSs, const Susy::Event &event,
                                           const DiLepEvtType &ll,
                                           const LeptonVector &ls)
{
  std::ostringstream oss;
  oss<<"MatrixPrediction passSrSs("<<(passSrSs?"true":"false")<<")"
     <<" "<<dilepDetails(event, ll, ls)
     <<" weight="<<m_weightComponents.fake;
  return oss.str();
}
//----------------------------------------------------------
std::string MatrixPrediction::hftTreeName() const
{
    std::string sn=sampleName();
    replace(sn, "period", "fake.period"); // avoid conflict with data trees
    return sn;
}
//----------------------------------------------------------
