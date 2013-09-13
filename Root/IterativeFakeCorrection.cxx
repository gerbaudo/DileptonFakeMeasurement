
#include "SusyTest0/IterativeFakeCorrection.h"

//-----------------------------------------------//
// Constructor
//-----------------------------------------------//
IterativeFakeCorrection::IterativeFakeCorrection() :
  myHist(),
  m_dbg(false),
  m_real(NULL),
  m_nIter(8)
{

  // Set the files
  string indir = "fakeplots/";
  //string append = "_Feb6_n0124_AltIso_rates";
  //string append = "_Feb12_n0124_rates";
  //string append = "_Feb18_n0127_altIso_rates";
  //string dappend = "_Feb18_n0127_rates";
  //string append = "_Feb17_n0117_rates";
  //string dappend = "_May2_n0139_rates";
  //string append = "_May2_n0139_rates";
  //string dappend = "_May6_n0139_NewElD0_NewMuIso_rates";
  //string append = "_May6_n0139_NewElD0_NewMuIso_rates";
  //string dappend = "_May24_n0139_rates";
  //string append = "_May24_n0139_rates";

  //string dappend = "_Jun25_n0139_MET_rates";
  //string append = "_Jun25_n0139_MET_rates";
  //string dappend = "_Jul8_n0144_rates";
  //string append = "_Jul8_n0144_NoMCMetCut_rates";
  //string dappend = "_Jul15_n0145_rates";
  //string append = "_Jul15_n0145_rates";
  
  //string dappend = "_Sep10_n0145_forDavide_rates";
  //string append = "_Sep10_n0145_forDavide_rates";  
  string dappend = "_Sep11_n0145_forDavide_rates";
  string append = "_Sep11_n0145_forDavide_rates";  

  // Data File
  m_data.file = new TFile((indir+"data"+dappend+".root").c_str());
  m_data.name = "Data";
  m_data.sname = "data";
  m_data.color = kBlack;
  m_data.marker = 20;

  // MC File
  m_mc.file = new TFile((indir+"mc"+append+".root").c_str());
  //m_mc.file = new TFile((indir+"heavyflavor"+append+".root").c_str());
  m_mc.name = "MC";
  m_mc.sname = "mc";
  m_mc.color = kRed;
  m_mc.marker = 25;

}
//-----------------------------------------------//
// Destructor
//-----------------------------------------------//
IterativeFakeCorrection::~IterativeFakeCorrection()
{

}

//-----------------------------------------------//
// Method to iterate
//-----------------------------------------------//
void IterativeFakeCorrection::iterate()
{

  // File
  //TFile* file = new TFile("corFake_altIso_2013.root","RECREATE");
  //TFile* file = new TFile("temp.root","RECREATE");
  //TFile* file = new TFile("corFake_altMuIso_May2_2013.root","RECREATE");
  //TFile* file = new TFile("corFake_May6_2013.root","RECREATE");
  //TFile* file = new TFile("corFake_May24_2013.root","RECREATE");
  //TFile* file = new TFile("corFake_Jun25_2013.root","RECREATE");
  //TFile* file = new TFile("corFake_Jul8_2013.root","RECREATE");
  //TFile* file = new TFile("corFake_Jul15_2013.root","RECREATE");
  //TFile* file = new TFile("corFake_Sep10_2013_forDavide.root","RECREATE");
  TFile* file = new TFile("corFake_Sep11_2013_forDavide.root","RECREATE");
  m_data.file->cd();

  // Place holders for histograms
  // [0] -- numerator -- tight
  // [1] -- denominator -- loose
  TH1F* data_low[2];
  TH1F* data_high[2];
  TH1F* mc_low[2];
  TH1F* mc_high[2];
  TH1F* corrected[2];

  vector<string> leps;
  leps.push_back("muon");
  leps.push_back("elec");

  for(uint il=0; il<leps.size(); ++il){
    string lep = leps.at(il);
    
    // Load Real Efficiency:
    TH1F* temp_num = getHist(m_data.file, lep+"_realCR_all_l_pt_num");
    TH1F* temp_den = getHist(m_data.file, lep+"_realCR_all_l_pt_den");
    if(m_dbg) cout << "Have Real: "<<temp_num<<" "<<temp_den<<endl;
    m_real = RatioHist(temp_num, temp_den, "");
    m_real->SetName("real_eff");
  
    // Load data
    data_low[0]  = getHist(m_data.file, lep+"_fakeHF_all_l_pt_num");
    data_low[1]  = getHist(m_data.file, lep+"_fakeHF_all_l_pt_den");
    data_high[0] = getHist(m_data.file, lep+"_fakeHF_high_all_l_pt_num");
    data_high[1] = getHist(m_data.file, lep+"_fakeHF_high_all_l_pt_den");

    // Load MC
    mc_low[0]  = getHist(m_mc.file, lep+"_fakeHF_all_l_pt_num");
    mc_low[1]  = getHist(m_mc.file, lep+"_fakeHF_all_l_pt_den");
    mc_high[0] = getHist(m_mc.file, lep+"_fakeHF_high_all_l_pt_num");
    mc_high[1] = getHist(m_mc.file, lep+"_fakeHF_high_all_l_pt_den");
    
    // Set the corrected. This is iteration 0,
    // so the constants are 0.
    corrected[0] = (TH1F*) data_low[0]->Clone("corrected_num");
    corrected[1] = (TH1F*) data_low[1]->Clone("corrected_den");
    
    // Now iterate and update the histogram
    for(int i=0; i<m_nIter; ++i){
      cout<<"--------------------------------"<<endl;
      cout<<"Iteration: "<<i<<endl;
      
      // Dump the rate
      for(int bin=1; bin<=corrected[0]->GetNbinsX(); ++bin){
	float num = corrected[0]->GetBinContent(bin);
	float den = corrected[1]->GetBinContent(bin);
	cout<<"Bin: "<<bin<<" num: "<<num<<" den: "<<den<<" rate: "<<num/den<<endl;
      }
      
      // Make temporary rate
      TH1F* rate = RatioHist(corrected[0], corrected[1], "");
      
      // Get Corrections
      vector<double> C_T = getC(rate, data_high[0], data_high[1], mc_high[0], true);
      vector<double> C_L = getC(rate, data_high[0], data_high[1], mc_high[1], false);
      
      // Correct the rates
      correctRate(corrected[0], data_low[0], mc_low[0], C_T);
      correctRate(corrected[1], data_low[1], mc_low[1], C_L);
      
      delete rate;
      
    }// end iteration
    
    
    TH1F* ratio = RatioHist(corrected[0],corrected[1],"");
    ratio->SetName((lep+"_corHFRate").c_str());
    file->cd();
    ratio->Write();
    m_data.file->cd();
    ratio->Delete();
    
    // Clean up
    data_low[0]->Delete();    data_low[1]->Delete();
    data_high[0]->Delete();   data_high[1]->Delete();
    mc_low[0]->Delete();      mc_low[1]->Delete();
    mc_high[0]->Delete();     mc_high[1]->Delete();
    corrected[0]->Delete();   corrected[1]->Delete();

  }// end loop over leptons

  // Write and close file
  file->Write();
  file->Close();



}

//-----------------------------------------------//
// Correct rate
//-----------------------------------------------//
void IterativeFakeCorrection::correctRate(TH1F* &rate, 
					  TH1F* data, 
					  TH1F* mc, 
					  vector<double> corrections)
{

  // Loop and correct each bin.  
  int nbins = rate->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    float d_bc = data->GetBinContent(bin);
    float d_be = data->GetBinError(bin);
    float m_bc = mc->GetBinContent(bin);
    float m_be = mc->GetBinError(bin);
    float corr = corrections.at(bin-1);
    cout<<"\t\tbin: "<<bin<<" data: "<<d_bc<<" mc: "<<m_bc<<" Corrected: "<<d_bc - corr * m_bc<<endl;
    rate->SetBinContent(bin, d_bc - corr * m_bc);
    rate->SetBinError(bin, sqrt(pow(d_be,2) + pow(corr*m_be,2)) );

  }

}

//-----------------------------------------------//
// Get correction
//-----------------------------------------------//
vector<double> IterativeFakeCorrection::getC(TH1F* rate,
					     TH1F* data_num,
					     TH1F* data_den,
					     TH1F* mc,
					     bool tight)
{

  // Correction factor:
  // C = (N^(data,high) - N^(fake pred, high))/N^(MC,high)

  cout<<"Getting corrections for tight? "<<tight<<endl;
  vector<double> corrections;

  int nbins = data_num->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){

    // Set N loose and tight for 1-D MM
    float nLoose = data_den->GetBinContent(bin);
    float nTight = data_num->GetBinContent(bin); //tight ? data_num->GetBinContent(bin) : 0.; 

    // Get r and f
    float r = m_real->GetBinContent(bin);
    float f = rate->GetBinContent(bin);
    //f = r-f <= 0.3 ? r-0.3 : f;

    // Set the values for C
    double nData = tight ? nTight : nLoose;
    double nMC   = mc->GetBinContent(bin);
    double nFake = getFake(r, f, nLoose, nTight, tight);
    cout<<"\tCorrection: "<<nData<<" "<<nFake<<" "<<nMC<<" === "<<(nData-nFake)/nMC<<endl;
    corrections.push_back( (nData - nFake)/nMC );

  }// end loop over bins

  return corrections;

}

//-----------------------------------------------//
// Get the fake contribution
//-----------------------------------------------//
double IterativeFakeCorrection::getFake(float r,
					float f,
					float nL,
					float nT,
					bool tight)
{

  cout<<"\t\t\tr: "<<r<<" f: "<<f<<" nL: "<<nL<<" nT: "<<nT<<endl;
  if(tight)
    return f/(r-f) * (nL*r - nT);
    
  return 1/(r-f) * (nL*r -nT);
}

//-----------------------------------------------//
TH1F* IterativeFakeCorrection::getHist(TFile* file, 
				       string histname)
{

  return (TH1F*) file->Get(histname.c_str());

}
