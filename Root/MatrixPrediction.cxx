#include <iomanip>
#include "SusyNtuple/SusyDefs.h"
#include "SusyTest0/MatrixPrediction.h"

using namespace std;
using namespace Susy;


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
  SusyPlotter()
{
}
//----------------------------------------------------------
void MatrixPrediction::Begin(TTree* /*tree*/)
{
  m_doFake = true;

  SusyPlotter::Begin(0);
  if(m_dbg) cout << "MatrixPrediction::Begin" << endl;
  // Load the matrix method package
  m_matrix = new SusyMatrixMethod::DiLeptonMatrixMethod();
  string pathRateFile = (string( std::getenv("ROOTCOREDIR"))
                         +"/../SusyMatrixMethod/data/pass6_Apr2_2013.root");
  m_matrix->configure(pathRateFile, SusyMatrixMethod::PT);
  cout<<"Matrix method initialized: "<<endl;
  bookFakeHisto();
  dump.open("fakeDump.txt");
}
//----------------------------------------------------------
Bool_t MatrixPrediction::Process(Long64_t entry)
{
  GetEntry(entry);
  clearObjects();
  if(m_do1fb && !is1fb()) return kTRUE;
  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%50000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << nt.evt()->run
         << " event " << setw(7) << nt.evt()->event << " ****" << endl;
  }
  selectObjects(); // select signal objects
  if( !selectEvent() )              return kTRUE;
  if( m_baseLeptons.size() != 2 )   return kTRUE;
  if( !passTrigger(m_baseLeptons) ) return kTRUE;
  bool count(true);
  SusyMatrixMethod::FAKE_REGION reg = SusyMatrixMethod::FR_VRSSbtag;
  SusyMatrixMethod::SYSTEMATIC  sys = SusyMatrixMethod::SYS_NONE;
  const Met*          m = m_met;
  const JetVector&    j = m_signalJets2Lep;
  const LeptonVector& l = m_baseLeptons;
  float metRel = getMetRel(m, l, j);
  float weight = getFakeWeight(l, reg, metRel, sys);
  // function references to shorten lines
//   void (&fh)(const LeptonVector &l, const JetVector &j, const Met* m,
//              const float weight, PlotRegion PR, uint sys) = fillHistos;
//   void (&ffh)(const LeptonVector &l, const JetVector &j, const Met *m,
//               float weight, PlotRegion PR, uint sys) = fillFakeHistos
  if(passSR6(l, j, m, count)) { fillHistos(l, j, m, weight, PR_SR6); fillFakeHistos(l, j, m, weight, PR_SR6, sys); }
  if(passSR7(l, j, m, count)) { fillHistos(l, j, m, weight, PR_SR7); fillFakeHistos(l, j, m, weight, PR_SR7, sys); }
  if(passSR8(l, j, m, count)) {
    cout<<"passSR8"<<endl;
    fillHistos(l, j, m, weight, PR_SR8); fillFakeHistos(l, j, m, weight, PR_SR8, sys); }
  if(passSR9(l, j, m, count)) {
    cout<<"passSR9"<<endl;
fillHistos(l, j, m, weight, PR_SR9); fillFakeHistos(l, j, m, weight, PR_SR9, sys); }
  return kTRUE;
}
//----------------------------------------------------------
void MatrixPrediction::Terminate()
{
  SusyPlotter::Terminate();
  if(m_dbg) cout << "MatrixPrediction::Terminate" << endl;
  delete m_matrix;
}
//----------------------------------------------------------
void MatrixPrediction::bookFakeHisto()
{
  m_histFile->cd();  // Histogram file from SusyPlotter
  for(uint iPR=0; iPR<PR_N; ++iPR){ // Plot Region names
    string PR = PRNames[iPR];
    //if( !(iPR == PR_VR1 || iPR == PR_VR3) ) continue;
    for(uint iCh=0; iCh<Ch_N; ++iCh){ // lepton channel loop
      string chan = chanNames[iCh];
      for(uint iMP=0; iMP<MP_N; ++iMP){
        string MP = MPNames[iMP];
        for(uint iWT=0; iWT<WTog_N; ++iWT){
          string WT = WTNames[iWT];
          for(uint iSYS=0; iSYS<m_systs.size(); ++iSYS){
            string SYS = m_systNames.at(iSYS);
            string base = PR+"_"+chan+"_"+MP+"_"+WT+"_"+SYS;
            // Preprocessor convenience: make a histogram by name (leave off the "h_") and binning
            #define NEWHIST(name, xLbl, nbin, min, max)			\
	      do{							\
		hf_ ## name[iCh][iPR][iMP][iWT][iSYS] = new TH1F((base+"_"+#name).c_str(), #name ";" xLbl, nbin, min, max); \
		hf_ ## name[iCh][iPR][iMP][iWT][iSYS]->Sumw2();		\
	      }while(0)

            #define NEWHIST2(name, xLbl, nbin, min, max)			\
	      do{							\
		hf_ ## name[iCh][iPR][iMP][iWT][iSYS] = new TH2F((base+"_"+#name).c_str(), #name ";" xLbl, nbin, min, max,nbin,min,max); \
		hf_ ## name[iCh][iPR][iMP][iWT][iSYS]->Sumw2();		\
	      }while(0)

            #define NEWVARHIST(name, xLbl, nbin, bins)				\
	      do{							\
		hf_ ## name[iCh][iPR][iMP][iWT][iSYS] = new TH1F((base+"_"+#name).c_str(), #name ";" xLbl, nbin, bins); \
		hf_ ## name[iCh][iPR][iMP][iWT][iSYS]->Sumw2();		\
	      }while(0)

            // Lepton Kin
            NEWHIST(l0_pt, "Lepton P_{T}", nptbins, ptmin, ptmax);
            NEWHIST(l1_pt, "Lepton P_{T}", nptbins, ptmin, ptmax);
            // Mass
            NEWHIST(ll_M, "m(ll)", nmassbins, massmin, massmax);
            // Met
            NEWHIST(met, "#slash{E}_{T}", nptbins, ptmin, ptmax);
            NEWHIST(metrel, "#slash{E}^{rel}_{T}", nptbins, ptmin, ptmax);
            // One bin
            NEWHIST(onebin, "onebin", 1, -0.5, 0.5);
            // event weight
            NEWHIST(evt_weight, "evt_weight", 5000, -5, 5);
            // l0 Pt vs l1 Pt
            NEWHIST2(pt0vspt1, "pt0vspt1", nptbins, ptmin, ptmax);
            // Mt plots
            NEWHIST(met_l0_Mt, "met_l0_Mt", nmassbins, massmin, massmax);
            NEWHIST(met_l1_Mt, "met_l1_Mt", nmassbins, massmin, massmax);
            // Jet plots
            NEWHIST(njets, "njets", 5, -0.5, 4.5);
            NEWHIST(bjet_pt, "bjet_pt", nptbins, ptmin, ptmax);
            NEWHIST(ljet_pt, "ljet_pt", nptbins, ptmin, ptmax);
           #undef NEWHIST
           #undef NEWHIST2
           #undef NEWVARHIST
          }// end loop over systemaitcs
        }// end loop over weight toggle
      }// end loop over Matrix Pairs
    }// end loop over Channel
  }// end loop over plot regions
}
//----------------------------------------------------------
void MatrixPrediction::fillFakeHistos(const LeptonVector &baseLeps, const JetVector &jets,
				      const Met* met,float weight, PlotRegion PR, uint sys)
{

  if(m_dbg) cout << "MatrixPrediction::plotFakeHisto" << endl;
  if( baseLeps.size() != 2 ) return;
  int ch = getChan(baseLeps);
  int mp = getMatrixPair(baseLeps);
  const Lepton* l0 = baseLeps[0];
  const Lepton* l1 = baseLeps[1];

  #define FILL(h, var)					\
    do{							\
      float max   = h[ch][PR][mp][WT_ON][sys]->GetXaxis()->GetXmax();	\
      float xfill = var > max ? max - 1e-4 : var;			\
      h[ch][PR][mp][WT_ON][sys]->Fill(xfill,weight);			\
      h[Ch_all][PR][mp][WT_ON][sys]->Fill(xfill,weight);		\
      h[ch][PR][mp][WT_OFF][sys]->Fill(xfill,1.0);			\
      h[Ch_all][PR][mp][WT_OFF][sys]->Fill(xfill,1.0);			\
      h[ch][PR][MP_ALL][WT_ON][sys]->Fill(xfill,weight);		\
      h[Ch_all][PR][MP_ALL][WT_ON][sys]->Fill(xfill,weight);		\
      h[ch][PR][MP_ALL][WT_OFF][sys]->Fill(xfill,1.0);			\
      h[Ch_all][PR][MP_ALL][WT_OFF][sys]->Fill(xfill,1.0);		\
    }while(0)

  #define FILL2(h, varx, vary)				\
    do{									\
      float maxx   = h[ch][PR][mp][WT_ON][sys]->GetXaxis()->GetXmax();	\
      float xfill  = varx > maxx ? maxx - 1e-4 : varx;			\
      float maxy   = h[ch][PR][mp][WT_ON][sys]->GetYaxis()->GetXmax();	\
      float yfill  = vary > maxy ? maxy - 1e-4 : vary;			\
      h[ch][PR][mp][WT_ON][sys]->Fill(xfill,yfill,weight);			\
      h[Ch_all][PR][mp][WT_ON][sys]->Fill(xfill,yfill,weight);		\
      h[ch][PR][mp][WT_OFF][sys]->Fill(xfill,yfill,1.0);			\
      h[Ch_all][PR][mp][WT_OFF][sys]->Fill(xfill,yfill,1.0);			\
      h[ch][PR][MP_ALL][WT_ON][sys]->Fill(xfill,yfill,weight);		\
      h[Ch_all][PR][MP_ALL][WT_ON][sys]->Fill(xfill,yfill,weight);		\
      h[ch][PR][MP_ALL][WT_OFF][sys]->Fill(xfill,yfill,1.0);			\
      h[Ch_all][PR][MP_ALL][WT_OFF][sys]->Fill(xfill,yfill,1.0);		\
    }while(0)

  float metrel = getMetRel(met, baseLeps, jets);
  FILL( hf_l0_pt, l0->Pt() );
  FILL( hf_l1_pt, l1->Pt() );
  FILL( hf_ll_M, (*l0 + *l1).M() );
  FILL( hf_met, met->Et );
  FILL( hf_metrel, metrel);
  FILL( hf_onebin, 0. );
  FILL( hf_evt_weight, weight);
  FILL2( hf_pt0vspt1, l0->Pt(), l1->Pt() );
  FILL( hf_njets, jets.size() );
  for(uint ij=0; ij<jets.size(); ++ij){
    Jet* jet = jets.at(ij);
    if( isCentralBJet(jet) ) FILL(hf_bjet_pt, jet->Pt());
    if( isCentralLightJet(jet) ) FILL(hf_ljet_pt, jet->Pt());
  }
  #undef FILL
  #undef FILL2
}
//----------------------------------------------------------
float MatrixPrediction::getFakeWeight(const LeptonVector &baseLeps,
                                      SusyMatrixMethod::FAKE_REGION region,
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
				      SusyMatrixMethod::FAKE_REGION region,
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
