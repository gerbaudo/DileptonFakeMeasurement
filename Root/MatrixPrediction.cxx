#include <iomanip>
#include "SusyNtuple/SusyDefs.h"
#include "SusyAna2012/MatrixPrediction.h"

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


/*--------------------------------------------------------------------------------*/
// Constructor
/*--------------------------------------------------------------------------------*/
MatrixPrediction::MatrixPrediction() : 
  SusyPlotter()
{

}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
/*--------------------------------------------------------------------------------*/
void MatrixPrediction::Begin(TTree* /*tree*/)
{
  m_doFake = true;

  SusyPlotter::Begin(0);
  if(m_dbg) cout << "MatrixPrediction::Begin" << endl;

  // Load the matrix method package
  m_matrix = new SusyMatrixMethod::DiLeptonMatrixMethod();
  cout<<"dir: "<<std::getenv("ROOTCOREDIR")<<endl;
  //string fakefile = "pass0_Moriond_Feb14_2013";
  //string fakefile = "pass1_Moriond_Feb15_2013";
  //string fakefile = "pass2_Moriond_Feb22_2013";
  //string fakefile = "pass3_Mar3_2013";
  //string fakefile = "test_Mar4_SherpaW";
  //string fakefile = "HtPtParam";
  //string fakefile = "Njets";
  //string fakefile = "test_VR3";
  //string fakefile = "test_AllWJets";
  //string fakefile = "test_NoQCDSF";
  //string fakefile = "test_FullSep_NoQCDSF";
  //string fakefile = "test_Mar5_NoQCDSF";
  //string fakefile = "test_Mar5_FullSep";
  //string fakefile = "test_Mar5_FullSep_SherpaW_AltMuIso";
  //string fakefile = "test_ChSep";
  //string fakefile = "pass4_Mar28_2013";
  //string fakefile = "pass5_Apr1_2013";
  //string fakefile = "pass6_Apr2_2013";
  //string fakefile = "test_NoSF_Apr2";
  //string fakefile = "test_Apr3";
  //string fakefile = "test_Apr3_altIso16";
  //string fakefile = "test_Apr4_MET";
  //string fakefile = "test_Apr5";
  //string fakefile = "test_Apr9_bjetParam_SherpaW";
  //string fakefile = "test_Apr9_SherpaW";
  //string fakefile = "test_Apr9_withBBbar";
  //string fakefile = "test_Apr9_L25";
  //string fakefile = "test_Apr9_withBBbar_L25_MV185";
  //string fakefile = "test_Apr9_altMuIso";
  //string fakefile = "test_Apr10_Sherpa";
  //string fakefile = "test_Apr30";
  //string fakefile = "test_Apr30_SherpaW";
  //string fakefile = "test_Apr30_d0";
  //string fakefile = "test_May6_newMuIso";
  //string fakefile = "test_May6_NewElD0_NewMuIso";
  //string fakefile = "test_May7";
  //string fakefile = "test_May7_ptParam";
  //string fakefile = "test_May8_MuEtconeCut";

  string fakefile = "pass0_Summer2013_May24";
  //string fakefile = "test_d0z0Baseline";
  //string fakefile = "test_d0z0tightBaseline";

  string cd = string(std::getenv("ROOTCOREDIR"));
  //string path = cd + "/../SusyMatrixMethod/data/"+fakefile+".root";
  string path = cd + "/../SusyPlotting/run/finalfake/"+fakefile+".root";
  cout<<"Path: "<<path<<endl;
  m_matrix->configure(path, 
		      SusyMatrixMethod::PT,     // Electron Real
		      SusyMatrixMethod::PT,     // Electron Fake
		      //SusyMatrixMethod::PT_ETA, // Electron Fake
		      SusyMatrixMethod::PT,     // Muon Real
		      SusyMatrixMethod::PT      // Muon Fake
		      );
  cout<<"Matrix method initialized: "<<endl;

  // book fake histograms
  bookFakeHisto();

  // file for dumping information
  dump.open("fakeDump.txt");

}

/*--------------------------------------------------------------------------------*/
// Main process loop function 
/*--------------------------------------------------------------------------------*/
Bool_t MatrixPrediction::Process(Long64_t entry)
{
  // Communicate tree entry number to SusyNtObject
  GetEntry(entry);
  clearObjects();

  if(m_do1fb && !is1fb()) return kTRUE;

  // Chain entry not the same as tree entry
  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%50000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << nt.evt()->run
         << " event " << setw(7) << nt.evt()->event << " ****" << endl;
  }

  // select signal objects
  selectObjects();

  // Check Analysis level cuts
  if( !selectEvent() )              return kTRUE;
  if( m_baseLeptons.size() != 2 )   return kTRUE;
  if( !passTrigger(m_baseLeptons) ) return kTRUE;

  if(m_doMCFake && !isFakeDilepton(m_baseLeptons)) return kTRUE;

  // Add additional checks to baseline leptons
  // Note this modifies fake rates!! make sure to use right ones
  //if( !isBaselineLepton(m_baseLeptons) ) return kTRUE; 


  //for(int s = 0; s<SusyMatrixMethod::SYS_N; ++s){
  //for(int s = 0; s<1; ++s){
  //SusyMatrixMethod::SYSTEMATIC sys = (SusyMatrixMethod::SYSTEMATIC) s;
  for(uint s = 0; s<m_systs.size(); ++s){
    SusyMatrixMethod::SYSTEMATIC sys = (SusyMatrixMethod::SYSTEMATIC) m_systs.at(s);

    // Plot Regions and use appropriate weights
    float metRel = getMetRel(m_met, m_baseLeptons, m_signalJets);

    //float mll = Mll(m_baseLeptons[0], m_baseLeptons[1]);
    //if( oppositeSign(m_baseLeptons) && fabs(mll - 91) < 15 ){
    //float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SROSjveto, metRel,sys);
    //fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_OSInc, s);
    //}
    /*
    // SS inclusive
    if( sameSign(m_baseLeptons) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SSInc, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SSInc,s);    
    }

    // OS Inc
    if( oppositeSign(m_baseLeptons) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SROSjveto, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_OSInc, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_OSInc,s);    
    }

    // SR1
    if( passSR1(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SROSjveto, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR1, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR1,s);    
    }
    // SR2
    if( passSR2(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRSSjets, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR2, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR2,s);
    }
    // SR3
    if( passSR3(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SR2jets, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR3, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR3,s);
    }
    // SR4
    if( passSR4(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRmt2a, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR4, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR4,s);
    }
    // SR4b
    if( passSR4b(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRmt2b, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR4b, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SR4b,s);
    }
    // SRWWa
    if( passWWSRa(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRWWa, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRWWa, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRWWa,s);
    }
    // SRWWb
    if( passWWSRb(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRWWb, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRWWb, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRWWb,s);
    }
    // SRWWc
    if( passWWSRc(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRWWc, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRWWc, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRWWc,s);
    }
    */
    // VR1
    if( passVR1(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR1, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR1,s);
    }
    // VR2
    /*
    if( passVR2(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSS, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR2, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR2,s);
    }
    */
    // VR3
    if( passVR3(m_baseLeptons, m_signalJets2Lep, m_met) ){
      //float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSS, metRel,sys);
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR3, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR3,s);
    }
    /*
    // VR4
    if( passVR4(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSS, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR4, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_VR4,s);
    }

    // WWCR1
    if( passWWCR1(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRWW1, metRel, sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_WWCR1, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_WWCR1, s);
    }
    // WWCR2
    if( passWWCR2(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRWW2, metRel, sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_WWCR2, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_WWCR2, s);
    }

    // TopCR
    if( passTOPCR(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRTOP, metRel, sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCR, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCR, s);
    }

    // TopCRWWa
    if( passTOPCRWWa(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRTOPWWa, metRel, sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCR, s);
      //fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCRWWa, s);
    }

    // TopCRWWb
    if( passTOPCRWWb(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRTOPWWb, metRel, sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCRWWb, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCRWWb, s);
    }
    // TopCRWWc
    if( passTOPCRWWa(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRTOPWWc, metRel, sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCRWWc, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_TOPCRWWc, s);
    }
    // Bonus region: Use VR2 weights
    if( passBR1(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR1, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR1,s);
    }
    if( passBR2(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR2, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR2,s);
    }
    if( passBR3(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR3, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR3,s);
    }
    if( passBR4(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR4, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR4,s);
    }
    if( passBR5(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_VRSSbtag, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR5, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_BR5,s);
    }
    */
    /*
    // ZJetSR
    if( passZJetsSR(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_SRZjets, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRZjets, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_SRZjets,s);
    }
    // ZJetCR-OS jet veto
    if( passZJetCRJVeto(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRZXOSjveto, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsJVeto, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsJVeto,s);
    }
    // ZJetCR-2jets
    if( passZJetCR2jets(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRZX2JETS, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjets2jets, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjets2jets,s);
    }
    // ZJetCR-mt2a
    if( passZJetCRMt2a(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRZXMT2a, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsMt2a, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsMt2a,s);
    }
    // ZJetCR-mt2b
    if( passZJetCRMt2b(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRZXMT2b, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsMt2b, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsMt2b,s);
    }
    // ZJetCR-WW
    if( passZJetCRWW(m_baseLeptons, m_signalJets2Lep, m_met) ){
      float weight = getFakeWeight(m_baseLeptons, SusyMatrixMethod::FR_CRZXWW, metRel,sys);
      fillHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsWW, s);
      fillFakeHistos(m_baseLeptons, m_signalJets2Lep, m_met, weight, PR_CRZjetsWW,s);
    }
    */
  }// end loop over systematics

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during a query
/*--------------------------------------------------------------------------------*/
void MatrixPrediction::Terminate()
{
  SusyPlotter::Terminate();
  if(m_dbg) cout << "MatrixPrediction::Terminate" << endl;
  
  delete m_matrix;
  
}

/*--------------------------------------------------------------------------------*/
// Book Fake Histos
/*--------------------------------------------------------------------------------*/
void MatrixPrediction::bookFakeHisto()
{

  // Histogram file from SusyPlotter
  m_histFile->cd();
  

  // Plot Region names
  for(uint iPR=0; iPR<PR_N; ++iPR){
    string PR = PRNames[iPR];
    //if( !(iPR == PR_VR1 || iPR == PR_VR2 || iPR == PR_VR3 || PR_TOPCR) ) continue;    
    if( !(iPR == PR_VR1 || iPR == PR_VR3) ) continue; 

    // lepton channel loop
    for(uint iCh=0; iCh<Ch_N; ++iCh){
      string chan = chanNames[iCh];
      
      for(uint iMP=0; iMP<MP_N; ++iMP){
	string MP = MPNames[iMP];

	for(uint iWT=0; iWT<WTog_N; ++iWT){
	  string WT = WTNames[iWT];

	  //for(int iSYS=0; iSYS<nFakeSys; ++iSYS){
	  //string SYS = FAKESYSNames[iSYS];
	  //for(int iSYS=0; iSYS<SusyMatrixMethod::SYS_N; ++iSYS){
	  //string SYS = SusyMatrixMethod::systematic_names[iSYS];
	  for(uint iSYS=0; iSYS<m_systs.size(); ++iSYS){
	    string SYS = m_systNames.at(iSYS);
	    
	    string base = PR+"_"+chan+"_"+MP+"_"+WT+"_"+SYS;

	    // Preprocessor convenience
	    // make a histogram by name (leave off the "h_") and binning
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

/*--------------------------------------------------------------------------------*/
// Plot the fake histograms
/*--------------------------------------------------------------------------------*/
void MatrixPrediction::fillFakeHistos(const LeptonVector &baseLeps, const JetVector &jets, 
				      const Met* met,float weight, PlotRegion PR, uint sys)
{

  if(m_dbg) cout << "MatrixPrediction::plotFakeHisto" << endl;

    // Get Channel for leptons
  // ** Only dealing with exactly two leptons 
  if( baseLeps.size() != 2 ) return;
  int ch = getChan(baseLeps);
  int mp = getMatrixPair(baseLeps);

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
  // Some useful Definitions
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
  
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

  // Lepton Kin
  FILL( hf_l0_pt, l0->Pt() );
  FILL( hf_l1_pt, l1->Pt() );

  // Mass
  FILL( hf_ll_M, (*l0 + *l1).M() );

  // Met related
  FILL( hf_met, met->Et );
  float metrel = getMetRel(met, baseLeps, jets);
  FILL( hf_metrel, metrel);

  FILL( hf_onebin, 0. );

  FILL( hf_evt_weight, weight);

  FILL2( hf_pt0vspt1, l0->Pt(), l1->Pt() );

  FILL( hf_met_l0_Mt, MyMt(*l0, met->lv()));
  FILL( hf_met_l1_Mt, MyMt(*l1, met->lv()));

  FILL( hf_njets, jets.size() );
  for(uint ij=0; ij<jets.size(); ++ij){
    Jet* jet = jets.at(ij);
    if( isCentralBJet(jet) ) FILL(hf_bjet_pt, jet->Pt());
    if( isCentralLightJet(jet) ) FILL(hf_ljet_pt, jet->Pt());
  }

  #undef FILL
  #undef FILL2

}

/*--------------------------------------------------------------------------------*/
// Get Fake Event weight based on region
/*--------------------------------------------------------------------------------*/
float MatrixPrediction::getFakeWeight(const LeptonVector &baseLeps, 
				      SusyMatrixMethod::FAKE_REGION region, 
				      float metRel,
				      SusyMatrixMethod::SYSTEMATIC sys)
{

  if(baseLeps.size() != 2) return 0.0;

  uint nVtx = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  
  //m_matrix->setDileptonType(baseLeps[0]->isEle(), baseLeps[1]->isEle());
  float weight = m_matrix->getTotalFake( isSignalLepton(baseLeps[0],m_baseElectrons, m_baseMuons,nVtx,isMC),
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

  if(!m_doMCFake) return weight;
  else return weight * getEvtWeight(baseLeps,true,true);

}
/*--------------------------------------------------------------------------------*/
// Needed for validation
/*--------------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------------*/
// Get Dilepton pair type
/*--------------------------------------------------------------------------------*/
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
