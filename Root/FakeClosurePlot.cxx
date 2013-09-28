// This plotting script will grow and be used to format
// histograms that can be shown in talks

#include "SusyPlotting/FancyPlotting.h"

//---------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------//
FancyPlotting::FancyPlotting(/*FPRunOption opt*/) : 
  myHist(),
  m_opt(RO_ALL),
  m_dbg(0),
  m_addIntegral(false),
  m_MCColor(kRed),
  m_makeTable(false)
{

  // Option, for checking later
  //m_opt = opt;


}

//---------------------------------------------------------------------//
// Initialize all files and things
//---------------------------------------------------------------------//
void FancyPlotting::init(FPRunOption opt)
{

  // Option, for checking later
  m_opt = opt;

  // Potential file objects
  File data;
  File top;
  File ZX;
  File Ztautau;
  File WW;
  File WZ;
  File ZZ;
  File qcd; 
  File WJet;
  //File Wgam;

  // For testing
  File Zjet;
  File Zdib;

  string inDir = "anaplots/";

  //string append     = "_May9_n0139_TrigTest.AnaHists";
  //string dataappend     = "_May9_n0139_TrigTest.AnaHists";
  //string fakeappend = "_May8_n0139_TrigTest.FakeHists"; 

  //string append     = "_May24_n0139.AnaHists";
  //string dataappend     = "_May24_n0139.AnaHists";
  //string fakeappend = "_May24_n0139.FakeHists"; 
  //string fakeappend = "__May29_n0139_d0z0Baseline.FakeHists"; 

  //string append     = "_May29_n0139_vr1Met.AnaHists";
  //string dataappend     = "_May29_n0139_vr1Met.AnaHists";
  //string fakeappend = "__May29_n0139_d0z0Baseline_vr1Met.FakeHists"; 

  //string append     = "_Jun7_n0139.AnaHists";
  //string dataappend     = "_Jun7_n0139.AnaHists";
  //string fakeappend = "__Jun7_n0139.FakeHists"; 

  //string append     = "_Jun18_n0139.AnaHists";
  //string dataappend = "_Jun18_n0139.AnaHists";
  //string fakeappend = "_Jun18_n0139.FakeHists"; 

  //string append     = "_Jun25_n0139.AnaHists";
  //string dataappend = "_Jun25_n0139.AnaHists";
  //string fakeappend = "_Jun25_n0139_bbbar_SF.FakeHists"; 
  //string fakeappend = "_Jun26_n0139_bbbar_SF.FakeHists"; 

  //string append     = "_Jul8_n0144.AnaHists";
  //string dataappend = "_Jul8_n0144.AnaHists";
  //string fakeappend = "_Jul8_n0144.FakeHists"; 
  //string fakeappend = "_Jul12_n0144.FakeHists"; 

  string append     = "_Jul17_n0145.AnaHists";
  string dataappend = "_Jul17_n0145.AnaHists";
  string fakeappend = "_Jul17_n0145.FakeHists"; 


  // Data
  data.file  = new TFile((inDir+"data"+dataappend+".root").c_str());
  //data.file  = new TFile((inDir+"data_1fb_Aug27.AnaHists.root").c_str());
  data.name  = "Data";
  data.sname = "dt";
  data.color = kBlack;
  data.ismc  = false;
  data.isfake = false;
    
  // top
  top.file  = new TFile((inDir+"top"+append+".root").c_str());
  top.name  = "Top";
  top.sname = "top";
  top.color = kRed+1;
  top.ismc  = true;
  top.isfake = false;
  
  // Z+X
  ZX.file  = new TFile((inDir+"ZX"+append+".root").c_str());
  ZX.name  = "Z+X";
  ZX.sname = "ZX";
  ZX.color = kOrange-2;
  ZX.ismc  = true;
  ZX.isfake = false;

  // Diboson
  WW.file  = new TFile((inDir+"WW"+append+".root").c_str());
  WW.name  = "WW";
  WW.sname = "WW";
  WW.color = kAzure + 4;
  WW.ismc  = true;
  WW.isfake = false;
  
  // Z->tautau
  //Ztautau.file = new TFile((inDir+"Ztautau"+append+".root").c_str());
  //Ztautau.name = "Z#tau#tau";
  //Ztautau.sname = "Ztautau";
  //Ztautau.color = kSpring+1;
  //Ztautau.ismc = true;
  //Ztautau.isfake = false;
  //Ztautau.xsLumi = 13;

  // QCD -- Matrix prediction
  qcd.file  = new TFile((inDir+"data"+fakeappend+".root").c_str());
  //qcd.file  = new TFile((inDir+"totalMCFake_Feb19_n0127_mcFake.MCFakeHists.root").c_str());
  qcd.name  = "Fake";
  qcd.sname = "fake";
  qcd.color = kGray;
  //qcd.ismc  = false;
  //qcd.isfake = false; //false;
  qcd.ismc  = true;
  qcd.isfake = true; //false;
  
  // Z+jet
  Zjet.file  = new TFile((inDir+"Zjet"+append+".root").c_str());
  Zjet.name  = "Z+jets";
  Zjet.sname = "Zjet";
  Zjet.color = kOrange-2;
  Zjet.ismc  = true;
  Zjet.isfake = false;

  // W+jet
  /*
  WJet.file  = new TFile((inDir+"WJet"+append+".root").c_str());
  WJet.name  = "W+jets";
  WJet.sname = "WJet";
  WJet.color = kMagenta;
  WJet.ismc  = true;
  WJet.isfake = false;
  */

  // WZ
  WZ.file  = new TFile((inDir+"WZ"+append+".root").c_str());
  //WZ.file  = new TFile((inDir+"lllnu_WZ_fix.AnaHists.root").c_str());
  WZ.name  = "WZ";
  WZ.sname = "WZ";
  WZ.color = kSpring+1;
  WZ.ismc  = true;
  WZ.isfake = false;

  // ZZ
  ZZ.file  = new TFile((inDir+"ZZ"+append+".root").c_str());
  ZZ.name  = "ZZ";
  ZZ.sname = "ZZ";
  ZZ.color = kOrange+8;
  ZZ.ismc  = true;
  ZZ.isfake = false;

  // Wgam
  //Wgam.file  = new TFile((inDir+"Wgamma"+append+".root").c_str());
  //Wgam.name  = "W#gamma";
  //Wgam.sname = "Wgam";
  //Wgam.color = kViolet;
  //Wgam.ismc  = true;
  //Wgam.isfake = false;

  
  m_files.clear();

  /*
  m_files.push_back(data);
  m_files.push_back(top);
  m_files.push_back(WW);
  m_files.push_back(ZX);
  m_files.push_back(Ztautau);
  m_files.push_back(qcd);
  */

  // Normal:
  m_files.push_back(data);
  m_files.push_back(ZX);
  //m_files.push_back(Wgam);
  m_files.push_back(top);
  m_files.push_back(WW);
  m_files.push_back(qcd);

  // Fake closure:
  //m_files.push_back(qcd);
  //m_files.push_back(top);
  //m_files.push_back(WJet);
  //m_files.push_back(WW);
  //m_files.push_back(ZX);


  //m_files.push_back(Ztautau); // in Zjet
  //m_files.push_back(ZX);

  /*
  m_files.push_back(data);
  m_files.push_back(qcd);
  m_files.push_back(Zdib);
  m_files.push_back(top);
  m_files.push_back(WW);
  //m_files.push_back(Ztautau); // in Zjet
  //m_files.push_back(ZX);
  m_files.push_back(Zjet);
  */


  if(opt == RO_ALL){
    m_PRs.push_back(PR_VRSS);
  }
  else return;

  setPlots();

}

//---------------------------------------------------------------------//
// Destructor
//---------------------------------------------------------------------//
FancyPlotting::~FancyPlotting()
{

}

//---------------------------------------------------------------------//
// Initialize vectors
//---------------------------------------------------------------------//
void FancyPlotting::setPlots()
{

  // Set any and all plots that are needed
  // and the functions will loop over this 
  // list

  // Lepton kinematics
  m_plots.push_back( pair<string,string> ("l0_pt", "l_{0} P_{T} [GeV]") );  
  m_plots.push_back( pair<string,string> ("l1_pt", "l_{1} P_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("e_pt", "Electron P_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("m_pt", "Muon P_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("m_eta", "#mu #eta") );  
  //m_plots.push_back( pair<string,string> ("l0_eta", "l_{0} #eta") );  
  //m_plots.push_back( pair<string,string> ("l1_eta", "l_{1} #eta") );

  // Mass plots
  m_plots.push_back( pair<string,string> ("ll_M", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("l0_j_M", "m(l0j) [GeV]") );
  //m_plots.push_back( pair<string,string> ("l1_j_M", "m(l1j) [GeV]") );
  //m_plots.push_back( pair<string,string> ("llj_M", "m(llj) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_fine", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_finer", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_optimal1", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_optimal2", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_optimal3", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_dPhiReg", "m(ll) [GeV] |d#phi|<3") );
  //m_plots.push_back( pair<string,string> ("ll_M_pos", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("ll_M_neg", "m(ll) [GeV]") );
  //m_plots.push_back( pair<string,string> ("llj_M", "m(llj) [GeV]") );
  //m_plots.push_back( pair<string,string> ("llj_M_pos", "m(llj) [GeV]") );
  //m_plots.push_back( pair<string,string> ("llj_M_neg", "m(llj) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_ll_Mt", "mt(ll,met) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_ll_Mt_onej", "mt(ll,met) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_ll_Mt_twoj", "mt(ll,met) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_ll_Mt_ge3j", "mt(ll,met) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_ll_Mt_oneOrtwoj", "mt(ll,met) [GeV]") );
  //m_plots.push_back( pair<string,string> ("dPhi_met_j", "dPhi(met,j)") );  
  //m_plots.push_back( pair<string,string> ("met_l0_Mt", "mt(met,l0) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_l1_Mt", "mt(met,l1) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_l0_Mt_eLead", "mt(met,e) [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_l0_Mt_muLead", "mt(met,m) [GeV]") );

  //m_plots.push_back( pair<string,string> ("j0_pt", "j_{0} P_{T}") );  
  //m_plots.push_back( pair<string,string> ("j0_eta", "j_{0} #eta") );  

  // Misc
  //m_plots.push_back( pair<string,string> ("sumQ", "Sum of Charge") );
  //m_plots.push_back( pair<string,string> ("njets_mll_90_120", "# jets") );
  //m_plots.push_back( pair<string,string> ("njets_mll_90_120_pos", "# jets") );
  //m_plots.push_back( pair<string,string> ("njets_mll_90_120_neg", "# jets") );
  //m_plots.push_back( pair<string,string> ("nbjets_mll_90_120", "# jets") );
  //m_plots.push_back( pair<string,string> ("nfjets_mll_90_120", "# jets") );
  //m_plots.push_back( pair<string,string> ("dR_llmet_j", "dR(llmet,jet)") );  


  //m_plots.push_back( pair<string,string> ("dPhi_llmet_j", "dPhi(llmet,jet)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_l0", "dPhi(met,l0)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_l1", "dPhi(met,l1)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_ll", "dPhi(met,ll)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_j", "dPhi(met,j)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_j0", "dPhi(met,j0)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_j1", "dPhi(met,j1)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_met_j2", "dPhi(met,j2)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_ll_j", "dPhi(ll,j)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_l0_j", "dPhi(l0,j)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_l1_j", "dPhi(l1,j)") );  
  //m_plots.push_back( pair<string,string> ("dPhi_l0_l1", "dPhi(l0,l1)") );
  //m_plots.push_back( pair<string,string> ("dPhi_met_bjet", "dPhi(met,bjet)") );
  //m_plots.push_back( pair<string,string> ("dPhi_met_bjet0", "dPhi(met,bjet0)") );
  //m_plots.push_back( pair<string,string> ("dPhi_met_bjet1", "dPhi(met,bjet0)") );
  //m_plots.push_back( pair<string,string> ("dPhi_met_bjet2", "dPhi(met,bjet0)") );
  //m_plots.push_back( pair<string,string> ("dPhi_met_bjet3", "dPhi(met,bjet0)") );

  /*
  m_plots.push_back( pair<string,string> ("dPhi_woSig_llmet_j", "dPhi(llmet,jet)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_met_l0", "dPhi(met,l0)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_met_l1", "dPhi(met,l1)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_met_ll", "dPhi(met,ll)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_met_j", "dPhi(met,j)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_ll_j", "dPhi(ll,j)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_l0_j", "dPhi(l0,j)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_l1_j", "dPhi(l1,j)") );  
  m_plots.push_back( pair<string,string> ("dPhi_woSig_l0_l1", "dPhi(l0,l1)") );  
  */


  // Met
  m_plots.push_back( pair<string,string> ("met", "#slash{E}_{T} [GeV]") );
  m_plots.push_back( pair<string,string> ("metrel", "#slash{E}^{rel}_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_refEle", "#slash{E}^{refEle}_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_refMuo", "#slash{E}^{refMuo}_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_refJet", "#slash{E}^{refJet}_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_softJet", "#slash{E}^{softJet}_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_refGamma", "#slash{E}^{refGamma}_{T} [GeV]") );
  //m_plots.push_back( pair<string,string> ("met_refCell", "#slash{E}^{cellOut}_{T} [GeV]") );



  // # jets
  m_plots.push_back( pair<string,string> ("njets", "# jets") );
  m_plots.push_back( pair<string,string> ("nbjets", "# b jets") );
  //m_plots.push_back( pair<string,string> ("bjet_pt", "CentralBJet P_{T}") );
  //m_plots.push_back( pair<string,string> ("ljet_pt", "CentralLighJet P_{T}") );

  // impact params

  //m_plots.push_back( pair<string,string> ("l0_d0", "l0 d0") );
  //m_plots.push_back( pair<string,string> ("l0_z0", "l0 z0") );
  //m_plots.push_back( pair<string,string> ("l0_d0sig", "l0 d0sig") );
  //m_plots.push_back( pair<string,string> ("l0_z0sintheta", "l0 z0*sin(#theta)") );
  //m_plots.push_back( pair<string,string> ("l1_d0", "l1 d0") );
  //m_plots.push_back( pair<string,string> ("l1_z0", "l1 z0") );
  //m_plots.push_back( pair<string,string> ("l1_d0sig", "l1 d0sig") );
  //m_plots.push_back( pair<string,string> ("l1_z0sintheta", "l1 z0*sin(#theta)") );
  //m_plots.push_back( pair<string,string> ("e_d0", "e d0") );
  //m_plots.push_back( pair<string,string> ("e_z0", "e z0") );
  //m_plots.push_back( pair<string,string> ("e_d0sig", "e d0sig") );
  //m_plots.push_back( pair<string,string> ("e_z0sintheta", "e z0*sin(#theta)") );
  //m_plots.push_back( pair<string,string> ("m_d0", "#mu d0") );
  //m_plots.push_back( pair<string,string> ("m_z0", "#mu z0") );
  //m_plots.push_back( pair<string,string> ("m_d0sig", "#mu d0sig") );
  //m_plots.push_back( pair<string,string> ("m_z0sintheta", "#mu z0*sin(#theta)") );
  /*
  m_plots.push_back( pair<string,string> ("l0_d0_biased", "l0 d0") );
  m_plots.push_back( pair<string,string> ("l0_z0_biased", "l0 z0") );
  m_plots.push_back( pair<string,string> ("l0_d0sig_biased", "l0 d0sig") );
  m_plots.push_back( pair<string,string> ("l0_z0sintheta_biased", "l0 z0*sin(#theta)") );
  m_plots.push_back( pair<string,string> ("l1_d0_biased", "l1 d0") );
  m_plots.push_back( pair<string,string> ("l1_z0_biased", "l1 z0") );
  m_plots.push_back( pair<string,string> ("l1_d0sig_biased", "l1 d0sig") );
  m_plots.push_back( pair<string,string> ("l1_z0sintheta_biased", "l1 z0*sin(#theta)") );
  m_plots.push_back( pair<string,string> ("e_d0_biased", "e d0") );
  m_plots.push_back( pair<string,string> ("e_z0_biased", "e z0") );
  m_plots.push_back( pair<string,string> ("e_d0sig_biased", "e d0sig") );
  m_plots.push_back( pair<string,string> ("e_z0sintheta_biased", "e z0*sin(#theta)") );
  m_plots.push_back( pair<string,string> ("m_d0_biased", "#mu d0") );
  m_plots.push_back( pair<string,string> ("m_z0_biased", "#mu z0") );
  m_plots.push_back( pair<string,string> ("m_d0sig_biased", "#mu d0sig") );
  m_plots.push_back( pair<string,string> ("m_z0sintheta_biased", "#mu z0*sin(#theta)") );
  */
  //m_plots.push_back( pair<string,string> ("n_preElecs", "# pre electrons") );
  //m_plots.push_back( pair<string,string> ("n_preMuons", "# pre muons") );

  clear();

}

//---------------------------------------------------------------------//
// Loop to make histograms
//---------------------------------------------------------------------//
void FancyPlotting::DataMCAnaPlots()
{

  // Here Loop over the channels and the signal regions that are set 
  // via the run options.  Right now we loop over all and save all.
  
  //string savedir = "formattedplots/May6_n0139_NewMuIso/";
  //string savedir = "formattedplots/May6_n0139_IPTurnedOff/";
  //string savedir = "formattedplots/May7_n0139_NewElD0_NewMuIso_PtParamOnly/";
  //string savedir = "formattedplots/May7_n0139_LowPtCut/";
  //string savedir = "formattedplots/May7_n0139_MuEtCone/";
  //string savedir = "formattedplots/May8_n0139_TrigTest/";
  //string savedir = "formattedplots/May24_n0139/";
  //string savedir = "formattedplots/May29_n0139/";
  //string savedir = "formattedplots/May29_n0139_vr1Met/";
  //string savedir = "formattedplots/May31_n0139/";
  //string savedir = "formattedplots/Jun7_n0139/";
  //string savedir = "formattedplots/Jun18_n0139/";
  //string savedir = "formattedplots/Jun25_n0139/";
  //string savedir = "formattedplots/Jun26_n0139/";
  //string savedir = "formattedplots/Jul8_n0144/";
  //string savedir = "formattedplots/Jul12_n0144/";
  string savedir = "formattedplots/Jul17_n0145/";

  for(uint ipr=0; ipr<m_PRs.size(); ++ipr){
    PlotRegion pr = m_PRs.at(ipr);
    string region = PRNames[pr];

    if(m_makeTable){
      dumpTable(pr);
      continue;
    }

    for(int ich=0; ich<Ch_N; ++ich){
      //if(ich != Ch_mm) continue;
      if(ich == Ch_all) continue;
      //if(ich == Ch_em) continue;
      Chan ch        = (Chan) ich;
      string ch_name = chanNames[ch];
      
      for(uint ip=0; ip<m_plots.size(); ++ip){
	string var    = m_plots.at(ip).first;
	string xtitle = m_plots.at(ip).second;

	// Get all histograms
	buildHists(m_hists, m_sys, var, xtitle, ch, pr);

	// Build Errors
	m_errs.push_back(buildErrors(m_hists.back(), m_sys));
	m_errs.push_back(buildRatioErrors(m_hists.back(), m_errs.back()));

	// Make the Legend
	float xleg[] = {0.65,0.9};
	float yleg[] = {0.5,0.89};
	TLegend* leg = buildLegend(m_hists, m_errs[0], xleg, yleg);

	// Plot
	string save = savedir + region + "_" + ch_name + "_" + var + ".pdf";
	//string save = savedir + region + "_" + ch_name + "_" + var + ".eps";
	//string save = savedir + region + "_" + ch_name + "_" + var + ".png";
	plotAll(m_hists, m_errs, save, leg, ich, false); //true);
	//plotAll(m_hists, m_errs, save, leg, ich, true);
	clear();


      }// end loop over plots
    }// end loop over channels
  }// end loop over regions

}

//---------------------------------------------------------------------//
// Set Histograms 
//---------------------------------------------------------------------//
void FancyPlotting::buildHists(vector<TH1F*> &hists, vector<TH1F*> &sys, string var, 
			       string xtitle, Chan ch, PlotRegion PR)			       
{
  if(m_dbg) cout << "FancyPlotting::buildHists from PR = " << PR << endl;

  //hists.clear();
  //sys.clear();
  
  string plot = PRNames[PR] + u() + chanNames[ch] + u() + var;
  if(m_dbg) cout << "Getting... " << plot <<endl;
  
  // Get a summary histogram to act as the total SM
  TH1F* SM = Get(m_files.at(0).file, plot + "_NOM", m_MCColor);
  SM->Reset();
  SM->SetName("SM");

  // Systematic histograms
  TH1F* sys_up = Get(m_files.at(0).file, plot + "_NOM", kBlack);
  sys_up->Reset();
  sys_up->SetName("sys_up");
  
  TH1F* sys_dn = Get(m_files.at(0).file, plot + "_NOM", kBlack);
  sys_dn->Reset();
  sys_dn->SetName("sys_dn");
  
  sys.push_back(sys_up);
  sys.push_back(sys_dn);

  // Loop over files and build histograms
  //cout<<"Plot Region: "<<PRNames[PR]<<" "<<chanNames[ch]<<endl;
  for(uint f=0; f<m_files.size(); ++f){

    File F = m_files.at(f);
    string end = !F.isfake ? "_NOM" : "_NONE";
    //string end = !F.ismc ? "_NOM" : "_NONE";
    //string end = f!=0 ? "_NOM" : "_NONE";
    //cout<<"File: "<<F.file<<" "<<plot+end<<endl;
    //string end = !F.isfake ? "" : "_NONE";
    TH1F* h = Get(F.file, plot + end, F.color, xtitle.c_str(), "Entries");
    h->SetName(F.name.c_str());
    h->SetFillColor(F.color);
    hists.push_back(h);

    // Get Sys
    if(F.isfake)    addFakeSys(h, F.file, plot, sys);
    else if(F.ismc) addSysError(h, F.file, plot, sys);

    // Add to SM
    if( F.ismc ) SM->Add(h, 1.);
    //cout<<"\tFile: "<<F.name<<" Integral: "<<Form("%4.2f",h->Integral())<<endl;
  }
  //cout<<"\tSM: "<<Form("%4.2f",SM->Integral())<<endl;

  //for(int bin=1; bin<=SM->GetNbinsX(); ++bin)
  //SM->SetBinError(bin, 0.);

  hists.push_back(SM);

  
}
//---------------------------------------------------------------------//
// Build legend
//---------------------------------------------------------------------//
TLegend* FancyPlotting::buildLegend(vector<TH1F*> hists, TGraphAsymmErrors* errs,
				    float* x, float* y)
{
  
  // All histograms that are in the list get
  // put into the legend. Will have to add 
  // the TGraphAsymmErrors later...
  
  // It is assumed that data is first and SM is last in the
  // histogram vector
  if(m_addIntegral){ x[0] -= 0.10; x[1] -= 0.10; }
  TLegend* legend = makeLegend(x,y);

  uint nhists = hists.size();
  if(m_addIntegral){
    for(uint i=0; i<nhists; ++i){
      float integral = hists.at(i)->Integral(0,-1);
      const char* name = hists.at(i)->GetName();
      if(i==0)             legend->AddEntry(hists.at(i), Form("%s %4.2f",name,integral), "P");
      else if(i==nhists-1) legend->AddEntry(hists.at(i), Form("%s %4.2f",name,integral), "L");
      else                 legend->AddEntry(hists.at(i), Form("%s %4.2f",name,integral),"F");
    }
  }
  else{
    for(uint i=0; i<nhists; ++i){
      if(i==0)             legend->AddEntry(hists.at(i), hists.at(i)->GetName(), "P"); 
      else if(i==nhists-1) legend->AddEntry(hists.at(i), hists.at(i)->GetName(), "L");
      else                 legend->AddEntry(hists.at(i), hists.at(i)->GetName(), "F");
    }
  }
  
  legend->AddEntry(errs, "Uncertainty","F");

  // --- In future, add Systematics here --- //
  
  return legend;

}

//---------------------------------------------------------------------//
// Build Ratio
//---------------------------------------------------------------------//
TH1F* FancyPlotting::buildRatio(TH1F* data, TH1F* SM)
{
  if(m_dbg) cout << "FancyPlotting::buildRatio" << endl;

  TH1F* ratio = RatioHist(data, SM, "Data/SM");
  return ratio;

}

//---------------------------------------------------------------------//
// Build the THStack from a set of histograms
//---------------------------------------------------------------------//
THStack* FancyPlotting::buildStack(vector<TH1F*> hists)
{
  if(m_dbg) cout << "FancyPlotting::buildStack" << endl;

  uint begin = hists.size() - 2; // Don't include total MC in the stack
  uint end   = 1;                // Don't include data in the stack

  THStack* stack = new THStack("stack","stack");
  for(uint i=begin; i>=end; i--)
    stack->Add( hists.at(i) );
  
  return stack;

}
//---------------------------------------------------------------------//
// Create and add error bars
//---------------------------------------------------------------------//
TGraphAsymmErrors* FancyPlotting::buildErrors(TH1F* summary, vector<TH1F*> sys)
{

  // Pass the summary histogram so we can get both the nominal 
  // values and the statistical error from it. Also will add
  // functions later to incorporate the various sources of error.

  const int s = 1000;
  double x[s];
  double xerr[s];
  double y[s];
  double yp[s];
  double ym[s];
  
  int nbins = summary->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    int ind = bin-1;
    xerr[ind] = summary->GetBinWidth(bin)/2.;
    x[ind]    = summary->GetBinCenter(bin);
    y[ind]    = summary->GetBinContent(bin);
    yp[ind]   = sqrt(sys[0]->GetBinContent(bin));
    ym[ind]   = sqrt(sys[1]->GetBinContent(bin));
    if(y[ind] - ym[ind] < 0) ym[ind] = y[ind];

    summary->SetBinError(bin, 0.); // don't want to double count this stat err
  }

  TGraphAsymmErrors* errs = new TGraphAsymmErrors(nbins, x, y, xerr, xerr, ym, yp);
  errs->SetMarkerSize(0);
  errs->SetFillStyle(3004);
  errs->SetFillColor(kGray+3);
  errs->SetLineWidth(2);

  return errs;
  
}
//---------------------------------------------------------------------//
TGraphAsymmErrors* FancyPlotting::buildRatioErrors(TH1F* nominal, TGraphAsymmErrors* tg_errs)
{

  // Pass the summary histogram so we can get both the nominal 
  // values and the statistical error from it. Also will add
  // functions later to incorporate the various sources of error.

  const int s = 1000;
  double x[s];
  double xerr[s];
  double y[s];
  double yp[s];
  double ym[s];
  
  int nbins = nominal->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    int ind = bin-1;
    float bc  = nominal->GetBinContent(bin);
    xerr[ind] = nominal->GetBinWidth(bin)/2.;
    x[ind]    = nominal->GetBinCenter(bin);

    Double_t xpos; Double_t ypos;
    tg_errs->GetPoint((Int_t)ind,xpos,ypos);
    float errDn = ypos - tg_errs->GetErrorYlow(ind);
    float errUp = ypos + tg_errs->GetErrorYhigh(ind);
    y[ind]    = 1.0;
    yp[ind]   = (bc !=0) ? fabs(errUp/bc - 1) : 0;
    ym[ind]   = (bc !=0) ? fabs(errDn/bc - 1)  : 0;
    //cout<<"bc: "<<bc<<" errUp: "<<errUp<<" errDn: "<<errDn<<endl;
    //cout<<"Error up: "<<yp[ind]<<" Error down: "<<ym[ind]<<endl;

    //if(y[ind] - ym[ind] < 0) ym[ind] = y[ind];

  }

  TGraphAsymmErrors* errs = new TGraphAsymmErrors(nbins, x, y, xerr, xerr, ym, yp);
  errs->SetMarkerSize(0);
  errs->SetFillStyle(3004);
  errs->SetFillColor(kGray+3);

  return errs;
  
}
//---------------------------------------------------------------------//
void FancyPlotting::addFakeSys(TH1F* nominal, TFile* file, string plot, 
			       vector<TH1F*> &sys)
{

  // There are 8 total sys shifts that have the following names:
  // * {EL,MU}_RE_{UP,DOWN} = 4 shifts
  // * {EL,MU}_FR_{UP,DOWN} = 4 shifts
  // * Also statistical error
  
  // Load Fake systematic shifts
  vector<TH1F*> shifts;

  shifts.push_back((TH1F*) file->Get((plot + "_EL_RE_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_EL_RE_DOWN").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_RE_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_RE_DOWN").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_EL_FR_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_EL_FR_DOWN").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_FR_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_FR_DOWN").c_str()));

  int nbins = nominal->GetNbinsX();  
  for(int bin=0; bin<=nbins; ++bin){
    float stat = pow(nominal->GetBinError(bin),2);
    float err_up = stat; 
    float err_dn = stat;
    float bc = nominal->GetBinContent(bin);
    for(uint s=0; s<shifts.size(); ++s){
      float shift = shifts.at(s)->GetBinContent(bin);
      if(shift > bc) err_up += pow(shift-bc,2);
      else           err_dn += pow(shift-bc,2);
    }// end loop over sys
    
    sys[0]->SetBinContent(bin, sys[0]->GetBinContent(bin) + err_up);
    sys[1]->SetBinContent(bin, sys[1]->GetBinContent(bin) + err_dn);
  }// end loop over bins
  
  // Clean up systematic histograms
  for(uint i=0; i<shifts.size(); ++i) shifts.at(i)->Delete();
  shifts.clear();
  
}
//---------------------------------------------------------------------//
void FancyPlotting::getFakeSys(TH1F* nominal, TFile* file, string plot, 
			       float &sysup, float &sysdn)
{

  // There are 8 total sys shifts that have the following names:
  // * {EL,MU}_RE_{UP,DOWN} = 4 shifts
  // * {EL,MU}_FR_{UP,DOWN} = 4 shifts
  // * Also statistical error
  
  // Load Fake systematic shifts
  vector<TH1F*> shifts;
  shifts.push_back((TH1F*) file->Get((plot + "_EL_RE_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_EL_RE_DOWN").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_RE_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_RE_DOWN").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_EL_FR_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_EL_FR_DOWN").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_FR_UP").c_str()));
  shifts.push_back((TH1F*) file->Get((plot + "_MU_FR_DOWN").c_str()));

  //int low = nominal->FindBin(100);
  //int high = nominal->FindBin(108);


  int nbins = nominal->GetNbinsX();  
  for(int bin=0; bin<=nbins; ++bin){
    float bc = nominal->GetBinContent(bin);
    for(uint s=0; s<shifts.size(); ++s){
      float shift = shifts.at(s)->GetBinContent(bin);
      if(shift > bc) sysup += pow(shift-bc,2);
      else           sysdn += pow(shift-bc,2);
    }// end loop over sys
  }// end loop over bins
  /*
  int nbins = nominal->GetNbinsX();  
  for(int bin=low; bin<=high; ++bin){
    //for(int bin=0; bin<=nbins; ++bin){
    float bc = nominal->GetBinContent(bin);
    for(uint s=0; s<shifts.size(); ++s){
      float shift = shifts.at(s)->GetBinContent(bin);
      if(shift > bc) sysup += pow(shift-bc,2);
      else           sysdn += pow(shift-bc,2);
    }// end loop over sys
  }// end loop over bins
  */
  
  sysup = sqrt(sysup);
  sysdn = sqrt(sysdn);

  // Clean up systematic histograms
  for(uint i=0; i<shifts.size(); ++i) shifts.at(i)->Delete();
  shifts.clear();
  
}
//---------------------------------------------------------------------//
void FancyPlotting::addSysError(TH1F* nominal, TFile* file, string plot,
				vector<TH1F*> &sys)
{

  // Will add the mc systematics later.  Right now
  // I am only adding the statistical error
  int nbins = nominal->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    float stat = pow(nominal->GetBinError(bin),2);
    float err_up = stat; 
    float err_dn = stat;
    
    //
    // Add Systematic loop here
    //

    sys[0]->SetBinContent(bin, sys[0]->GetBinContent(bin) + err_up);
    sys[1]->SetBinContent(bin, sys[1]->GetBinContent(bin) + err_dn);
  }// end loop over bins    
    

}

//---------------------------------------------------------------------//
// General tools
//---------------------------------------------------------------------//
void FancyPlotting::plotAll(vector<TH1F*> hists, vector<TGraphAsymmErrors*> errs,
			    string save, TLegend* leg, int ch, bool logy, bool logx)			    
{

  if(m_dbg) cout << "plotAll" << endl;

  // TPad and Virtual Pad
  TCanvas* c = makeCanvas("c");
  TPad* _pTop;
  TPad* _pBot;
  makePads(c,_pTop,_pBot);
  _pTop->SetBottomMargin(0.02);
  _pBot->SetTopMargin(0.015);

  if(logx){
    _pTop->SetLogx();
    _pBot->SetLogx();
  }

  // Plotting three basic objects:
  TH1F* data     = hists.at(0);
  THStack* stack = buildStack(hists);
  TH1F* SM       = hists.at(hists.size()-1);
  TH1F* ratio    = buildRatio(data, SM);

  // Create lines for ratio plot
  TLine* nom = getLine(data, 1.00, kBlack, 1);
  TLine* up50 = getLine(data, 1.5, kBlack, 2);
  TLine* dn50 = getLine(data, 0.5, kBlack, 2);
  /*
  TLine* up1 = getLine(data, 2, kBlack, 2);
  TLine* up2 = getLine(data, 3, kBlack, 2);
  TLine* up3 = getLine(data, 4, kBlack, 2);
  TLine* up4 = getLine(data, 5, kBlack, 2);
  TLine* up5 = getLine(data, 6, kBlack, 2);
  */
  // Set the attributes for the top Pad
  // via the data histogram
  c->cd();
  _pTop->Draw();
  _pTop->cd();
  stack->Draw("hist");
  stack->GetYaxis()->SetTitleSize(0.065);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.06); //,"Y");
  stack->GetXaxis()->SetLabelSize(0.0); //, "X");
  stack->GetYaxis()->SetLabelFont(42);
  vector<float> minmax =  getMinMax(hists);

  if(!logy){
    stack->SetMinimum(minmax[0]);
    stack->SetMaximum(minmax[1]);
  }
  else if( logy ){
    _pTop->SetLogy();
    stack->SetMinimum(1e-4);
    stack->SetMaximum(minmax[1] * 100);
  }

  // Draw the top canvas objects
  stack->Draw("hist");
  SM->Draw("same hist");
  errs[0]->Draw("E2 same");
  data->Draw("p same pE");
  leg->Draw("same");

  TLatex* lat = makeLatex();
  lat->SetTextSize(0.055);
  lat->DrawLatex(0.4, 0.85, "#intLdt=21fb^{-1}");
  string chname = "";
  if( ch == Ch_ee) chname = "ee";
  if( ch == Ch_mm) chname = "#mu#mu";
  if( ch == Ch_em) chname = "e#mu";
  lat->DrawLatex(0.45, 0.75, chname.c_str());
  _pTop->Update();

  // Fix ratio plots to show +/- 100%
  ratio->SetMaximum(2.0);
  ratio->SetMinimum(0.0);
  ratio->SetLabelFont(42,"X");
  ratio->SetLabelFont(42,"Y");
  ratio->GetYaxis()->CenterTitle();
  ratio->GetXaxis()->SetTitleOffset(1.0);
  ratio->GetXaxis()->SetTitleSize(0.12);

  // Fix Errors. Current bounds 0 and 2
  //for(int bin=1; bin<=ratio->GetNbinsX(); ++bin){
  //float bc = ratio->GetBinContent(bin);
  //    float be = 
    
  // Draw the ratio plot on the bottom
  // Will need to add the error object
  c->cd();
  _pBot->Draw();
  _pBot->cd();
  ratio->Draw("ep");
  errs[1]->Draw("E2 same");
  nom->Draw("same");
  up50->Draw("same");
  dn50->Draw("same");
  //up1->Draw("same");
  //up2->Draw("same");
  //up3->Draw("same");
  //up4->Draw("same");
  //up5->Draw("same");
  ratio->Draw("ep same");
  // Save -- save is full path and name for plot
  c->SetFillColor(kWhite);
  c->Update();
  c->SaveAs( save.c_str() );

  delete c;
  delete lat;

}
//----------------------------------------------------//
float FancyPlotting::getNorm(TH1* h)
{
  return h->Integral(0,-1);
}
//-----------------------------------------------------//
float FancyPlotting::getStat(TH1* h, float low, float high)
{

  int nbins = h->GetNbinsX();
  float be = 0;
  int min = h->FindBin(low);
  int max = h->FindBin(high);
  for(int bin=min; bin<=max; ++bin){
    be += pow(h->GetBinError(bin),2);
  }
  
  return sqrt(be);

}
//-----------------------------------------------------//
void FancyPlotting::setMinMax(TH1* &h, float min, float max)
{
  h->SetMinimum(min);
  h->SetMaximum(max);
}
//-----------------------------------------------------//
TLine* FancyPlotting::getLine(TH1* h, float y, int color, int style)
{
  float x0     = h->GetBinCenter(1) - h->GetBinWidth(1)/2.;
  int finalBin = h->GetNbinsX();
  float x1     = h->GetBinCenter(finalBin) + h->GetBinWidth(finalBin)/2.;
  TLine* line  = makeLine(x0, x1, y, y, color, style);
  return line;
}
//-----------------------------------------------------//
float FancyPlotting::getMax(TH1F* h[], int n)
{
  float max = -999;
  for(int i=0; i<n; ++i){
    if( max < h[i]->GetMaximum() )
      max = h[i]->GetMaximum();
  }
  
  return max;

}
 
//-----------------------------------------------------//
// Add-on -- Dumping to a table
//-----------------------------------------------------//
void FancyPlotting::dumpTable(PlotRegion reg)
{

  // So this is kind of last minute, but need to be 
  // able to dump the validation region plots
  // We will need to get the fake error

  cout<<"Table for: "<<PRNames[reg]<<endl;

  
  // Header for table
  //------------------------------------//
  cout<<"Process & ee & (stat) & (sys) &";
  cout<<"$\\mu\\mu$ & (stat) & (sys) &";
  cout<<"e$\\mu$ & (stat) & (sys)\\\\";
  cout<<endl;
  //------------------------------------//

  vector<string> channels; 
  channels.push_back("ee");
  channels.push_back("mm");
  channels.push_back("em");

  float total[3]={0,0,0};
  float totstat[3]={0,0,0};
  float totsysup[3]={0,0,0};
  float totsysdn[3]={0,0,0};

  float low = 100;
  float high = 108;

  for(uint f=0; f<m_files.size(); ++f){
    File file = m_files.at(f);
    if( !file.ismc && !file.isfake ) cout<<"Observed ";
    else cout<<file.name;

    for(uint i=0; i<channels.size(); ++i){    
      cout<<" & ";

      //string plot = PRNames[reg] + "_" + channels[i] + "_ll_M_optimal3";
      string plot = PRNames[reg] + "_" + channels[i] + "_onebin";

      /*
      string plot = PRNames[reg] + "_" + channels[i] + "_onebin";
      if(reg == PR_VR1 && i == 0) // ee with Z veto
	plot = PRNames[PR_VR4] + "_" + channels[i] + "_onebin";
      if(reg == PR_VR2 && i == 0) // ee with Z veto
	plot = PRNames[PR_VR3] + "_" + channels[i] + "_onebin";
      */

      TH1F* nominal = NULL;
      if(file.isfake)
	nominal = (TH1F*) file.file->Get((plot+"_NONE").c_str());
      else
	nominal = (TH1F*) file.file->Get((plot+"_NOM").c_str());

      float nom  = nominal->GetBinContent(1);
      float stat = nominal->GetBinError(1);

      //float binl = nominal->FindBin(low);
      //float binh = nominal->FindBin(high);
      //float nom  = nominal->Integral(binl, binh);
      //float stat = getStat(nominal, low, high);

      float sysup = 0;
      float sysdn = 0;
      if(file.isfake) getFakeSys(nominal, file.file, plot, sysup, sysdn);

      if(file.ismc || file.isfake){
	total[i] += nom;
	totstat[i] += stat*stat;
	totsysup[i] += sysup*sysup;
	totsysdn[i] += sysdn*sysdn;
      }
      cout<<Form("%4.2f",nom)<<" & "
	  <<Form("%4.2f",stat)<<" & ";
      if(file.isfake) cout<<"$^{+"<<Form("%4.2f",sysup)<<"}_{-"<<Form("%4.2f",sysdn)<<"}$ ";
      else cout<<" $--$ ";
    }// end loop over channels
    cout<<" \\\\"<<endl;
    cout<<"\\hline"<<endl;
  }// end loop over files
  
  cout<<"Total SM ";
  for(uint i=0; i<channels.size(); ++i){
    cout<<" & "<<Form("%4.2f",total[i])
	<<" & "<<Form("%4.2f",sqrt(totstat[i]))
	<<" & $^{+"<<Form("%4.2f",sqrt(totsysup[i]))
	<<"}_{-"<<Form("%4.2f",sqrt(totsysdn[i]))<<"}$ ";
  }
  cout<<" \\\\"<<endl;
  cout<<"\\hline"<<endl;
		      


}
