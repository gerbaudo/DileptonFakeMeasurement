#include "SusyTest0/FakePlotting.h"
#include "SusyTest0/utils.h"

#include <cassert>
#include <iterator>

//--------------------------------------------------------//
// Constructor
//--------------------------------------------------------//
FakePlotting::FakePlotting(RunOption runOpt) :
  myHist(),
  m_dbg(0),
  m_runopt(RO_N),
  m_fakeStat(NULL)
{

  m_runopt = runOpt;
  string cd   = string(std::getenv("ROOTCOREDIR"));
  //string fakefile = "fakeRate_trial9_Nov2";
  //string fakefile = "pass1_Moriond_Feb15_2013";
  //string fakefile = "pass2_Moriond_Feb22_2013";
  string fakefile = "pass6_Apr2_2013";
  //string fakefile = "test_Feb25_FullSep";
  //string path = cd + "/../SusyMatrixMethod/data/" + fakefile + ".root";
  string path = cd + "/../SusyTest0/run/finalfake/" + fakefile + ".root";
  //m_fakeStat = new FakeStatTool::StatErrorTool(path,SusyMatrixMethod::PT);


}

//--------------------------------------------------------//
// Initialize all files
//--------------------------------------------------------//
bool FakePlotting::initInputFile(File &out,
                                 const string &fname, const string &name, const string &sname,
                                 int color, int marker)
{
  out.file  = new TFile(fname.c_str());
  out.name  = name;
  out.sname = sname;
  out.color = color;
  out.marker = marker;
  bool inputIsValid(out.file && out.file->IsOpen());
  if(!inputIsValid)
    cout<<"invalid input '"<<fname<<"'"<<" ("<<out.file<<")"<<endl;
  return inputIsValid;
}
//----------------------------------------------------------
FakePlotting& FakePlotting::setTag(const std::string &name)
{
  m_tag = name;
  return *this;
}
//----------------------------------------------------------
FakePlotting& FakePlotting::setInputDir(const std::string &dir)
{
  if(!dirExists(dir)) cout<<"Warning, invalid input dir '"<<dir<<"'"<<endl;
  m_inputdir = dir;
  return *this;
}
//----------------------------------------------------------
FakePlotting& FakePlotting::setOuputDir(const std::string &dir)
{
  m_outputdir = mkdirIfNeeded(dir);
  bool dirIsInvalid(m_outputdir.size()==0);
  if(dirIsInvalid) {
    string fallbackOutdir("./");
    cout<<"FakePlotting::setOuputDir('"<<dir<<"') : failed to create the output directory"<<endl
        <<"output files will be in '"<<fallbackOutdir<<"'"<<endl;
    m_outputdir = fallbackOutdir;
  }
  return *this;
}
//----------------------------------------------------------
FakePlotting& FakePlotting::setInputItercorrFile(const std::string &name)
{
  m_inputItercorrFile = name;
  return *this;
}
//----------------------------------------------------------
void FakePlotting::init()
{

  string inDir = m_inputdir+"/";
  initInputFile(data,    inDir+"data"       +m_tag+".root", "Data",        "data",  kBlack,    20);
  initInputFile(totMC,   inDir+"allBkg"     +m_tag+".root", "MC Combined", "mc",    kRed,      25);
  initInputFile(ttbar,   inDir+"ttbar"      +m_tag+".root", "t#bar{t}",    "ttbar", kRed,      21);
  initInputFile(Zjets,   inDir+"zjets"      +m_tag+".root", "Z+jets",      "Zjets", kOrange,   23);
  initInputFile(Wjets,   inDir+"wjets"      +m_tag+".root", "W+jets",      "Wjets", kViolet+2, 22);
  initInputFile(diboson, inDir+"diboson"    +m_tag+".root", "Diboson",     "dib",   kMagenta,  34);
  initInputFile(HF,      inDir+"heavyflavor"+m_tag+".root", "HF",          "HF",    kBlue,     26);

  if(m_runopt == RO_Data){
    m_files.push_back(data);
  }

  if(m_runopt == RO_MC || m_runopt == RO_SRComp || m_runopt == RO_DumpPer){
    cout<<"Pushing back files..."<<endl;
    m_files.push_back(HF);
    m_files.push_back(ttbar);
    m_files.push_back(Wjets);
    m_files.push_back(Zjets);
    m_files.push_back(diboson);
    totMC.color = kBlack;
    //if(m_runopt != RO_SRComp) m_files.push_back(totMC);
  }
  if(m_runopt == RO_DataMCSub || 
     m_runopt == RO_DataMCSF  ||
     m_runopt == RO_DataMC ){
    m_files.push_back(data);
    m_files.push_back(totMC);
  }
  if(m_runopt == RO_GJetCR){
    m_files.push_back(data);
    m_files.push_back(gjet);
  }

}

//--------------------------------------------------------//
// Destructor
//--------------------------------------------------------//
FakePlotting::~FakePlotting()
{

  // Clean up the files
  for(uint i=0; i<m_files.size(); ++i)
    m_files.at(i).file->Delete();
  m_files.clear();

  delete m_fakeStat;

}

//--------------------------------------------------------//
// Loop over and plot fake rate vs. various variables
//--------------------------------------------------------//
void FakePlotting::MCFakeRate()
{
  if(m_dbg) cout << "MCFakeRate" << endl;
  string savedir = m_outputdir;
  vector<string> plots, labels;
  plots.push_back("l_pt");   labels.push_back("P_{T}");
  TH1F* rates[6];
  TCanvas* c = makeCanvas("c");
  float x[] = {0.15, 0.3};
  float y[] = {0.90 - m_files.size()/20., 0.90};
  float Min_Max[] = {0,1.1};
  float fakeMinMax[] = {0,1.1};
  float realMinMax[] = {0.5, 1.3};
  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];
    for(uint ip=0; ip<plots.size(); ++ip){
      string var    = plots.at(ip);
      for(int cr =0; cr<CR_N; ++cr){
        if(!(cr == CR_MCConv || cr == CR_MCReal || cr == CR_MCQCD) ) continue;
        if( cr == CR_MCReal ){ Min_Max[0] = realMinMax[0]; Min_Max[1] = realMinMax[1];}
        else{ Min_Max[0] = fakeMinMax[0]; Min_Max[1] = fakeMinMax[1];}
        string crname = CRNames[cr];
        string xlabel = (il == LT_EL ? "Electron " : "Muon " )+labels.at(ip);
        string pname = lepname + "_" + crname + "_all_" + var;
        if(m_dbg) cout << "Grabbing: " << pname << endl;
        vector<string> h_names;
        for(uint f=0; f<m_files.size(); ++f){
          File F = m_files.at(f);
          rates[f] = buildRate(F.file, pname, F.sname, xlabel,"Rate", F.color, F.marker);
          h_names.push_back(F.name);
        }
        vector<float> frac = buildFraction(m_files, pname);
        cout<<"Pushing back the labels"<<endl;
        vector<Label> lbls;
        lbls.push_back( makeLabel(CRLabels[cr]     , 0.6, 0.87) );
        cout<<"building fractions"<<endl;
        cout<<"plotting"<<endl;
        string save = savedir + pname + ".pdf";
        plotFake(rates, h_names, c, Min_Max, x, y, lbls, save);
        for(uint f=0; f<m_files.size(); ++f) rates[f]->Delete();
      }// end loop over control region
    }// end loop over plots
  }// end loop over lepton types
}
//--------------------------------------------------------//
// Loop over and plot fake rate vs. various variables
//--------------------------------------------------------//
void FakePlotting::DataFakeRate()
{
  if(m_dbg) cout << "DataFakeRate" << endl;

  string savedir = mkdirIfNeeded("out/fakeplot/dataRate/");

  vector<string> plots;      vector<string> labels;
  plots.push_back("l_pt");   labels.push_back("P_{T}");
  //plots.push_back("l_eta");  labels.push_back("#eta");
  
  TH1F* rates[6];
  TCanvas* c = makeCanvas("c");
  
  float x[] = {0.15, 0.3};
  float y[] = {0.90 - CR_N/20., 0.90};

  float Min_Max[] = {0,1.1};

  vector<int> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kGreen);

  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];

    for(uint ip=0; ip<plots.size(); ++ip){
      string var    = plots.at(ip);
      string xlabel = labels.at(ip);
      
      vector<string> h_names;
      int nplots = 0;

      // Loop over the control regions and
      // plot the data rates from the various 
      // region son the same plot.
      for(int cr=0; cr<CR_N; ++cr){
	string crname = CRNames[cr];
	string crtitle = CRLabels[cr];

	string pname = lepname + "_" + crname + "_" + var;
	
	cout << "Grabbing: " << pname << endl;
	
    File F = m_files.at(0);
    TH1F *rate = buildRate(F.file, pname, F.sname, xlabel, "Rate", colors.at(nplots));
    if(rate) { rates[nplots] = rate; h_names.push_back(crtitle); nplots++; }
    else { cout<<"cannot build rate histo, skipping"<<endl; continue; }

	vector<Label> lbls;
	//lbls.push_back( makeLabel(CRLabels[cr],  0.5, 0.9) );
	
	string save = savedir + lepname + "_" + crname + "_" + var + ".eps";
	
	plotFake(rates, h_names, c, Min_Max, x, y, lbls, save);
      }// end loop over control region

    }// end loop over plots
  }// end loop over lepton types
  
}
//--------------------------------------------------------//
// Data and Photon rates
//--------------------------------------------------------//
void FakePlotting::GammaJetCRRates()
{
  if(m_dbg) cout << "GammaJetCRRates" << endl;

  string savedir = mkdirIfNeeded("out/fakeplot/gammajetCR/");

  vector<string> plots;      vector<string> labels;
  plots.push_back("l_pt");   labels.push_back("P_{T}");
  /*
  TH1F* rates[6];
  TCanvas* c = makeCanvas("c");
  
  float x[] = {0.15, 0.3};
  float y[] = {0.90 - m_files.size()/20., 0.90};

  float Min_Max[] = {0,1.1};

  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];

    for(uint ip=0; ip<plots.size(); ++ip){
      string var    = plots.at(ip);
      string xlabel = labels.at(ip);
      
      vector<string> h_names;
      string crname = CRNames[CR_FD];
      string pname = lepname + "_" + crname + "_" + var;
	
      cout << "Grabbing: " << pname << endl;
      
      for(uint f=0; f<m_files.size(); ++f){
	File F = m_files.at(f);
	rates[f] = buildRate(F.file, pname, F.sname, "Rate", xlabel, F.color);
	h_names.push_back(F.name);
      }
	
      vector<Label> lbls;
      lbls.push_back( makeLabel(CRLabels[CR_FD],  0.5, 0.87) );
      
      string save = savedir + pname + ".eps";
      
      plotFake(rates, h_names, c, Min_Max, x, y, lbls, save);
	

    }// end loop over plots
  }// end loop over lepton types
  */
}

//--------------------------------------------------------//
// Loop over and plot fake rate vs. various variables
//--------------------------------------------------------//
void FakePlotting::DataRateMCSub()
{
  cout<<"Needs fixed!!!!"<<endl;
  return;
  if(m_dbg) cout << "DataFakeRate" << endl;

  // Only need Data and MC files (so 2 total!)
  if( m_files.size() != 2 ){
    cout << "Not enough files to data - mc rates" << endl;
    return;
  }

  string savedir = mkdirIfNeeded("out/fakeplot/corDataRate/");
  
  File f_data = m_files[0];
  File f_mc   = m_files[1];
  
  vector<string> plots;      vector<string> labels;
  plots.push_back("l_pt");   labels.push_back("P_{T}");
  //plots.push_back("l_eta");  labels.push_back("#eta");
  
  TH1F* rates[2];
  TCanvas* c = makeCanvas("c");
  
  float x[] = {0.15, 0.3};
  float y[] = {0.90 - 2/20., 0.90};

  float Min_Max[] = {0,1.1};

  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];
    string lepton  = il == LT_EL ? "Electron" : "Muon";

    for(uint ip=0; ip<plots.size(); ++ip){
      string var    = plots.at(ip);
      string xlabel = lepton + " " + labels.at(ip);

      for(int cr=0; cr<CR_N; ++cr){

	//if(cr == CR_RA || cr == CR_FD || cr == CR_FE ) continue;

	string crname = CRNames[cr];
      
	string pname = lepname + "_" + crname + "_" + var;
	
	cout << "Grabbing: " << pname << endl;

	vector<string> h_names;	

	//rates[0] = buildRate(f_data.file, pname, f_data.sname, xlabel, "Rate", kBlack);
	//rates[1] = getCorrectedRate(f_data.file, f_mc.file, pname, xlabel, kRed);
	h_names.push_back("Data");
	h_names.push_back("Data - Contamination (MC)");
	
	vector<Label> lbls;
	lbls.push_back( makeLabel(CRLabels[cr],  0.5, 0.87) );	
	
	string save = savedir + pname + ".eps";
      
	plotFake(rates, h_names, c, Min_Max, x, y, lbls, save);
      }// end loop over control region	

    }// end loop over plots
  }// end loop over lepton types

}
//--------------------------------------------------------//
// Get SF for Data and MC
//--------------------------------------------------------//
void FakePlotting::DataMCSF(RunOption ro /*DG: 'ro' now unused? try to drop it*/)
{
  if(m_dbg) cout << "DataFakeRate" << endl;
  if( m_files.size() != 2 ){ // Only need Data and MC files (so 2 total!)
    cout << "Not enough files to data/mc sf" << endl;
    return;
  }
  string savedir = m_outputdir;
  vector<string> plots, labels;
  plots.push_back("all_l_pt"        ); labels.push_back("P_{T}"     );
  plots.push_back("all_l_eta"       ); labels.push_back("|#eta|"    );
  plots.push_back("all_l_eta_coarse"); labels.push_back("|#eta|"    );
  plots.push_back("all_onebin"      ); labels.push_back("Single Bin");
  TH1F* rates[3];
  TCanvas* c = makeCanvas("c");
  float x[] = {0.15, 0.3};
  float y[] = {0.90 - 2/20., 0.90};
  float Min_Max[] = {0,1.1};
  vector<Label> lbls;
  vector<string> names;
  TFile* f_temp = new TFile(m_inputItercorrFile.c_str());
  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];
    string lepton  = il == LT_EL ? "Electron" : "Muon";
    for(uint ip=0; ip<plots.size(); ++ip){
      string var    = plots.at(ip);
      string xlabel = lepton + " " + labels.at(ip);
      string ylabel = "Fake Rate";
      // For this I decided to opt out of a loop, since each CR needs different things.
      lbls.clear();
      names.clear();
      // Z tag and probe
      if(m_dbg) cout << "\tTrying Real T and P" << endl;
      string realZ = lepname + "_realCR_" + var;
      rates[0] = buildSideBandSubRate(data.file, lepname, var, xlabel, "Real eff", data.color, data.marker);
      rates[1] = buildSideBandSubRate(totMC.file, lepname, var, xlabel, "Real eff", totMC.color, totMC.marker);
      names.push_back("Data: Z Tag and Probe");
      names.push_back("MC Comb: Z Tag and Probe");
      plotSF(rates, names, c, Min_Max, x, y, lbls, savedir + "sf_" + realZ + ".pdf");
      names.clear();
      rates[0]->Delete();
      rates[1]->Delete();
      // HF tag and probe; in this case I need to subract out the contamination
      if(m_dbg) cout << "\tTrying HF" << endl;
      string tpHF     = lepname + "_fakeHF_" + var;
      string trueHF   = lepname + "_fakeHF_" + var + "_heavy";
      string contamHF   = lepname + "_fakeHF_" + var + "_others";
      cout<<"Getting: "<<tpHF<<" and for mc: "<<trueHF<<endl;
      // In this case totMC doesn't contain bb/cc
      rates[0] = (TH1F*) f_temp->Get((lepname+"_corHFRate").c_str());
      rates[0]->SetLineColor(data.color);
      rates[0]->SetMarkerColor(data.color);
      rates[0]->SetMarkerStyle(data.marker);
      rates[1] = buildRate(HF.file, tpHF, HF.sname, xlabel, ylabel, HF.color, HF.marker);
      names.push_back("Data: HF Tag and Probe (Iterative Subtraction)");
      names.push_back("b#bar{b}/c#bar{c} MC: HF Tag and Probe");
      cout<<"Going to plot HF sf: "<<endl;
      dumpHisto(rates[0]);
      dumpHisto(rates[1]);
      plotSF(rates, names, c, Min_Max, x, y, lbls, savedir + "sf_" + tpHF + ".pdf");
      names.clear();
      rates[0]->Delete();
      rates[1]->Delete();
      if(lepname == "elec"){ // Conversion
        if(m_dbg) cout << "\tTrying Conv" << endl;
        string conv = lepname + "_fakeConv_" + var;
        string convTrue = lepname + "_fakeConv_" + var;
        rates[0] = buildRate(data.file, conv, data.sname, xlabel, ylabel, data.color, data.marker);
        rates[1] = buildRate(totMC.file, conv, totMC.sname, xlabel, ylabel, totMC.color, totMC.marker);
        names.push_back("Data: Conv CR");
        names.push_back("MC Comb: Conv CR");
        plotSF(rates, names, c, Min_Max, x, y, lbls, savedir + "sf_" + conv + ".pdf");
        names.clear();
        rates[0]->Delete();
        rates[1]->Delete();
      }
    }// end loop over plots
  }// end loop over lepton types
  f_temp->Close();
}
//--------------------------------------------------------//
// Get true MC rates in signal regions
//--------------------------------------------------------//
void FakePlotting::MCSRRates()
{

  // This will plot the true fake rate from each signal region
  // so that we can compre to see if it is worth splitting up
  // fake rate into individual control regions.
  // Later maybe we can use it to find out what contributes
  // the most to each control region

  if(m_dbg) cout << "MCSRRates" << endl;

  string savedir = mkdirIfNeeded("out/fakeplot/mcSRRates/");

  vector<string> plots;      vector<string> labels;
  plots.push_back("l_pt");   labels.push_back("P_{T}");
  //plots.push_back("l_eta");  labels.push_back("#eta");
  
  TH1F* rates[6];
  TCanvas* c = makeCanvas("c");
  
  float x[] = {0.15, 0.3};
  float y[] = {0.90 - SR_N/20., 0.90};

  float Min_Max[] = {0,1.1};

  vector<int> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kViolet+2);
  colors.push_back(kOrange);
  colors.push_back(kMagenta);

  vector<string> h_names;

  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];

    for(int ls=0; ls<LS_N; ++ls){
      string lsname = LSNames[ls];
    //for(int ls=0; ls<LSS_N; ++ls){
    //string lsname = LSSNames[ls];

      for(uint ip=0; ip<plots.size(); ++ip){
	string var    = plots.at(ip);
	string xlabel = (il == LT_EL ? "Electron " : "Muon " )+labels.at(ip);

	h_names.clear();	

	for(int cr=0; cr<SR_N-1; ++cr){
	  string crname = SRNames[cr];
	  
	  string pname = lepname + "_" + crname + "_" + lsname + "_" + var;
	  
	  if(m_dbg) cout << "Grabbing: " << pname << endl;	  

	  
	  rates[cr] = buildRate(totMC.file, pname, crname, xlabel, "Rate", colors.at(cr), 20 + cr);
	  h_names.push_back(crname);

	}// end loop over signal regions
	
	vector<Label> lbls;

	string save = savedir + lepname + "_" + lsname + "_" + var + ".eps";
	
	plotFake(rates, h_names, c, Min_Max, x, y, lbls, save);

      }// end loop over plots
    }// end loop over LSSNames
  }// end loop over lepton types
  

}

//--------------------------------------------------------//
// Make the plot of composition in the signal regions
//--------------------------------------------------------//
void FakePlotting::Composition()
{
  
  // Plot the composition (HF/LF/Conversion) in the various
  // signal regions.
  
  string savedir = mkdirIfNeeded("out/fakeplot/Composition2012_Sep25/");

  TH1F* hists[6];
  TCanvas* c = makeCanvas("c");
  
  float x[] = {0.60, 0.90};
  float y[] = {0.95 - m_files.size()/20., 0.95};

  float Min_Max[] = {0,1.1};

  for(int il=0; il<LT_N; ++il){
    string lepname = LTNames[il];

    for(int sr=0; sr<SR_N; ++sr){
      string srname = SRNames[sr];
	  
      string pname = lepname + "_flavorComp_" + srname;
      
      if(m_dbg) cout << "Grabbing: " << pname << endl;
      
      vector<string> h_names;
      
      for(uint f=0; f<m_files.size(); ++f){
	File F = m_files.at(f);
	hists[f] = Get(F.file, pname, F.color,"","",F.marker);
	hists[f]->SetFillColor(F.color);
	//hists[f]->GetXaxis()->SetBinLabel(1,"HF");
	//hists[f]->GetXaxis()->SetBinLabel(2,"LF");
	//hists[f]->GetXaxis()->SetBinLabel(3,"Conv");
	//rates[f] = buildRate(F.file, pname, F.sname, xlabel, "Rate", F.color, F.marker);
	h_names.push_back(F.name);
      }
      
      vector<Label> lbls;
      lbls.push_back( makeLabel(SRNames[sr], 0.2, 0.87) );

      string save = savedir + pname + ".eps";
      
      vector<float> min_max = getMinMax(hists, h_names.size());
      Min_Max[0] = min_max[0];
      Min_Max[1] = min_max[1];
      plotStack(hists, h_names, c, x, y, lbls, save, true);
      
    }// end loop over signal region
  }// end loop over lepton types
  

}

//--------------------------------------------------------//
// dump SR composition
//--------------------------------------------------------//
void FakePlotting::dumpSRTable()
{

  for(int lt=0; lt<LT_N; ++lt){
  //for(int lt=0; lt<1; ++lt){
    string lep = LTNames[lt];
    cout<<"--------------------------------------"<<endl;
    cout<<"Lepton:  "<<lep<<endl;

    //cout<<"Signal Region & Heavy [\\%] & Light [\\%] & Conversion [\\%] \\\\"<<endl;
    cout<<"Signal Region & QCD [\\%] & Conversion [\\%] \\\\"<<endl;
    cout<<"\\hline"<<endl;
    for(int sr=0; sr<SR_N; ++sr){
      //for(int sr=0; sr<4; ++sr){
      string srname = SRNames[sr];
      // DG 2013-09-16 this is really an ugly hack! to be removed ASAP.
      if(SR_WHSS==sr) srname = CRNames[CR_SRWHSS];

      string pname = lep + "_" + srname + "_all_flavor_den";
      //string pname = lep + "_VR2_flavor_den";
    
      if(m_dbg) cout << "Grabbing: " << pname << endl;
      
      TH1F* temp = Get(totMC.file, pname);
      if(!temp) {
        cout<<"dumpSRTable: missing histo '"<<pname<<"' ...skip"<<endl;
        continue;
      }
      //TH1F* temp = Get(Wjets.file, pname);
      
      //float qcd   = temp->GetBinContent(1) + temp->GetBinContent(2);
      float qcd   = temp->GetBinContent(LS_QCD+1);
      //float HF    = temp->GetBinContent(LS_HF+1);
      //float LF    = temp->GetBinContent(LS_LF+1);
      float conv  = temp->GetBinContent(LS_Conv+1);
      //float total = HF + LF + conv;
      float total = qcd + conv;

      //cout<<"HF: "<<HF<<" LF: "<<LF<<" Conv: "<<conv<<endl;

      cout<<SRProperNames[sr]<<" & ";
      //cout<<Form("%4.0f",(HF/total *100))<<" & ";
      //cout<<Form("%4.0f",(LF/total *100))<<" & ";
      cout<<Form("%4.0f",(qcd/total *100))<<" & ";
      cout<<Form("%4.0f",(conv/total *100))<<"\\\\"<<endl;
      cout<<"\\hline"<<endl;
    }
  }



}

//--------------------------------------------------------//
// Make plots showing tight and loose info
//--------------------------------------------------------//
void FakePlotting::TLPlots()
{  
  // was using fakePred.file, Matt says it is not needed
}

//--------------------------------------------------------//
// Dump Signal Region prediction for fake
//--------------------------------------------------------//
void FakePlotting::dumpSRFake()
{
  // was using fakePred.file. Matt says it is not needed
}

//--------------------------------------------------------//
// Test method to getting additional stat error
//--------------------------------------------------------//
float FakePlotting::getNewStat(TFile* file, string sr, string ch)
{

  string base = sr + "_" + ch;
  float current = getBE(file, base + "_ALL_wOn_NONE_onebin");


  FakeStatTool::DileptonType chan;
  if( ch == "ee" )      chan = FakeStatTool::DT_ee;
  else if( ch == "em" ) chan = FakeStatTool::DT_em;
  else if( ch == "mm" ) chan = FakeStatTool::DT_mm;
  else return current;

  SusyMatrixMethod::FAKE_REGION reg;
  if( sr == "srmt2a" ) reg = SusyMatrixMethod::FR_SRmT2a;
  else if( sr == "srmt2b" ) reg = SusyMatrixMethod::FR_SRmT2b;
  else if( sr == "srmt2c" ) reg = SusyMatrixMethod::FR_SRmT2c;
  else if( sr == "srwwa" ) reg = SusyMatrixMethod::FR_SRWWa;
  else if( sr == "srwwb" ) reg = SusyMatrixMethod::FR_SRWWb;
  else if( sr == "srwwc" ) reg = SusyMatrixMethod::FR_SRWWc;
  else if( sr == "srzjets" ) reg = SusyMatrixMethod::FR_SRZjets;
  else return current;

  int nTT = getBC(file, base + "_TT_wOff_NONE_onebin");
  int nTL = getBC(file, base + "_TL_wOff_NONE_onebin");
  int nLT = getBC(file, base + "_LT_wOff_NONE_onebin");
  int nLL = getBC(file, base + "_LL_wOff_NONE_onebin");

  float newError = m_fakeStat->getUpdatedError(current,
					       nTT,
					       nTL,
					       nLT,
					       nLL,
					       reg,
					       chan);

  //cout<<"\t\t\tCombos: "<<nTT<<" "<<nTL<<" "<<nLT<<" "<<nLL<<endl;
  //cout<<"\t\t\tOld: "<<current<<" New: "<<newError<<endl;

  return newError;
}
				     
  

//--------------------------------------------------------//
// Attempt to get some HF Normalization...
//--------------------------------------------------------//
void FakePlotting::GetHFNorm()
{

  // Get the contributions from all of the various 
  // MC sources.

  vector<string> leptons;
  leptons.push_back("elec");
  leptons.push_back("muon");

  for(uint l=0; l<leptons.size(); ++l){

    string plot = leptons.at(l) + "_fakeHF_l_pt_den";
    string p_data = leptons.at(l) + "_fakeHF_l_pt_den";
    
    float n_wjet   = ((TH1F*) Wjets.file->Get(plot.c_str()))->Integral();
    float n_zjet   = ((TH1F*) Zjets.file->Get(plot.c_str()))->Integral();
    float n_ttbar  = ((TH1F*) ttbar.file->Get(plot.c_str()))->Integral();
    float n_dib    = ((TH1F*) diboson.file->Get(plot.c_str()))->Integral();
    float n_heavy  = ((TH1F*) HF.file->Get(plot.c_str()))->Integral();
    float n_data = ((TH1F*) data.file->Get(p_data.c_str()))->Integral();
    cout<<"--------------------"<<endl;
    cout<<leptons.at(l)<<endl;
    cout<<"W+jet:      "<<n_wjet<<endl;	
    cout<<"Z+jet:      "<<n_zjet<<endl;
    cout<<"ttbar:      "<<n_ttbar<<endl;
    cout<<"diboson:    "<<n_dib<<endl;
    cout<<"HF:         "<<n_heavy<<endl;
    cout<<"data:       "<<n_data<<endl;

    float norm = (n_data-n_wjet-n_zjet-n_ttbar-n_dib)/n_heavy;
    cout<<"Norm:       "<<norm<<endl;
  }
    
}

//--------------------------------------------------------//
// Plotting Functions
//--------------------------------------------------------//
void FakePlotting::plotFake(TH1F* h[], vector<string> names, TCanvas* c, 
			    float* MinMax, float* xLeg, float* yLeg,
			    vector<Label> lbls, string save, bool restrictMax, 
			    bool logy)
{
  if(m_dbg) cout << "plotFake with " << names.size() << " histograms"  << endl;
  c->cd();

  // Make Legend and add entries
  TLegend* legend = makeLegend(xLeg, yLeg);
  legend->SetTextSize(0.03);

  for(uint i=0; i<names.size(); ++i)
    legend->AddEntry(h[i],names.at(i).c_str(),"lep");

  // Make Latex object for the labels
  TLatex* lat = makeLatex();
  lat->SetTextSize(0.03);
  

  // Set Canvas Attributes
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.12);

  // Set Histogram attributes
  float max = getMaximum(h, names.size()) * 1.2;
  float min = getMinimum(h, names.size());
  if( min < 0 ) min = min; 
  else min *= 0.6;

  min = MinMax[0];
  max = MinMax[1];

  //min = min < 0 ? 0 : min;
  if(restrictMax) h[0]->SetMaximum( max > 1.3 ? 1.3 : max );
  else if(logy); // do nothing
  else h[0]->SetMaximum(max); 
  if(!logy) h[0]->SetMinimum(min);
  h[0]->GetYaxis()->SetTitleOffset(0.9);
  h[0]->GetYaxis()->SetTitleSize(0.05);
  h[0]->SetLabelSize(0.04, "X");
  h[0]->SetLabelSize(0.04, "Y");
  

  h[0]->Draw("ep");
  for(uint i=0; i<names.size(); ++i)
    h[i]->Draw("same ep");
  legend->Draw("same");
  for(uint l=0; l<lbls.size(); ++l)
    lat->DrawLatex( lbls[l].x, lbls[l].y, lbls[l].lbl.c_str() );
  
  //if(logy) c->SetLogy();
  c->SaveAs(save.c_str());

}
//--------------------------------------------------------//
// Plot SF
//--------------------------------------------------------//
void FakePlotting::plotSF(TH1F* h[], vector<string> names, TCanvas* c,
                          float* MinMax, float* xLeg, float* yLeg,
                          vector<Label> lbls, string save, bool doFit,
                          bool topBot, string topLabel, string ratLabel)
{
  // It is assumed here that h[0] is data and  h[1] is mc and the ratio will be Data/MC
  // Consider changing this later so we can have more plots on one canvas
  if( names.size() > 4 ){
    cout<<"Error in plotSF"<<endl<<"nHisto = "<<names.size()<<endl<<"Require <= 4, returning"<<endl;
    return;
  }
  TLegend* leg = makeLegend(xLeg,yLeg);
  for(uint i = 0; i<names.size(); ++i) leg->AddEntry(h[i],names.at(i).c_str(),"P");
  TPad *pTop, *pBot;
  TVirtualPad* tv = CreatePad(c,pTop,pBot,0.4);
  TLatex* lat = makeLatex();
  // Take care of drawing top first
  pTop->Draw();
  pTop->cd();
  h[0]->GetYaxis()->SetTitleSize(0.065);
  h[0]->GetYaxis()->SetTitleOffset(0.85);
  h[0]->SetLabelSize(0.06,"Y");
  vector<float> min_max = getMinMax(h, names.size());
  h[0]->SetMaximum(min_max[1] * 1.4);
  h[0]->SetMinimum(min_max[0] * 0.6);
  h[0]->GetYaxis()->SetTitle( topLabel.c_str() );
  h[0]->Draw("ep");
  for(uint i=1; i<names.size(); ++i) h[i]->Draw("same ep");
  for(uint i=0; i<lbls.size(); ++i){ Label l = lbls.at(i); lat->DrawLatex(l.x,l.y,l.lbl.c_str()); }
  leg->Draw("same");
  pTop->Update();
  // Bottom objects are the ratio with respect to the first plot. It is assumed there is atleast 2
  TH1F* ratio[4];
  for(uint i =0; i<names.size()-1; ++i){
    int color = h[i+1]->GetLineColor();
    assert(topBot); // othewise obsolete (DG 2013-09-26)
    ratio[i] = RatioHist(h[0],  h[1],ratLabel.c_str(),color);
  }
  ratio[0]->GetXaxis()->SetTitleSize(0.12);
  ratio[0]->GetXaxis()->SetTitleOffset(0.8);
  ratio[0]->SetLabelSize(0.09, "X");
  ratio[0]->GetYaxis()->SetTitleSize(0.08);
  ratio[0]->GetYaxis()->SetTitleOffset(0.7);
  ratio[0]->SetLabelSize(0.09,"Y");
  ratio[0]->GetYaxis()->CenterTitle();
  ratio[0]->GetYaxis()->SetNdivisions(510);

  TH1F* temp[1]; temp[0] = ratio[0]; // set min max
  min_max = getMinMax(temp,1);
  float MIN = getMinimum(ratio, names.size()-1);
  float MAX = getMaximum(ratio, names.size()-1);
  ratio[0]->SetMinimum( (1-MIN) < 0.1 ? 0.9 : MIN * 0.5 );
  ratio[0]->SetMaximum( (MAX-1) < 0.1 ? 1.1 : MAX * 1.2 );
  // Take care of drawing bottom
  tv->cd();
  pBot->Draw();
  pBot->cd();
  // Fit the ratio by a constant function
  float fitmin = ratio[0]->GetXaxis()->GetXmin();
  float fitmax = ratio[0]->GetXaxis()->GetXmax();
  TF1 f = TF1("func","[0]",fitmin,fitmax);
  f.SetParameter(0,1.);
  f.SetLineWidth(0);
  vector<float> fitval;
  vector<float> fiterr;
  vector<float> chi2;
  vector<int> ndf;
  if(doFit){
    for(uint i=0; i<names.size()-1; ++i){
      ratio[i]->Fit("func","RQ");
      ratio[i]->Draw();
      fitval.push_back(f.GetParameter(0));
      fiterr.push_back(f.GetParError(0));
      chi2.push_back(f.GetChisquare());
      ndf.push_back(f.GetNDF());
    }
  }
  ratio[0]->Draw("ep");
  for(uint i=1; i<names.size(); ++i) ratio[i-1]->Draw("same ep");
  if(doFit){
    TLatex* l = makeLatex();
    for(uint i = 0; i<names.size()-1; ++i){
      int color = ratio[i]->GetLineColor();
      l->DrawLatex(0.15 + i*0.25, 0.40, Form("#color[%i]{Fit =  %4.4f +/- %4.4f}",color,fitval.at(i), fiterr.at(i)));
      l->DrawLatex(0.15 + i*0.25, 0.32, Form("#color[%i]{#chi^{2}/#dof = %4.2f/%i}",color,chi2.at(i),ndf.at(i)));
    }
  }
  if(isCanvasWithNeededSf(save)) cout<<"SF for '"<<save<<"' : "<<vfloat2str(fitval)<<endl;
  c->SaveAs(save.c_str());
  replace(save, ".pdf", ".png");
  c->SaveAs(save.c_str());
  for(uint i=0; i<names.size()-1; ++i) ratio[i]->Delete();
}
//--------------------------------------------------------//
// Plot Stack
//--------------------------------------------------------//
void FakePlotting::plotStack(TH1F* h[], vector<string> names, TCanvas* c, 
			     float* xLeg, float* yLeg,
			     vector<Label> lbls, string save, bool logy)
{
  c->cd();

  // Make the stack
  THStack* stack = buildStack(h, names.size());
  
  // Make the legend
  TLegend* leg = makeLegend(xLeg,yLeg);
  leg->SetTextSize(0.05);
  float qcd  = 0;
  float conv = 0;
  int nentries = 0;
  for(uint i = 0; i<names.size(); ++i){
    qcd += h[i]->GetBinContent(1) + h[i]->GetBinContent(2);
    conv += h[i]->GetBinContent(4);
    nentries += h[i]->GetEntries();
    leg->AddEntry(h[i],names.at(i).c_str(),"F");
  }

  float total = qcd+conv;
  cout<<"Total entries: "<<nentries<<" %qcd: "<<qcd/total<<" %conv: "<<conv/total<<endl;  

  // Draw and set margins
  stack->Draw("hist");
  stack->GetYaxis()->SetTitleSize(0.065);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.06); //,"Y");
  stack->GetXaxis()->SetLabelSize(0.06); //, "X");
  stack->GetYaxis()->SetLabelFont(42);
  //for(int bin = 1; bin<= h[0]->GetNbinsX(); ++bin){
  //cout<< h[0]->GetXaxis()->GetBinLabel(bin) << endl;
  //stack->GetXaxis()->SetBinLabel(bin, h[0]->GetXaxis()->GetBinLabel(bin));
  //}  

  vector<float> minmax =  getMinMax(h, names.size());

  if(!logy){
    stack->SetMinimum(minmax[0]);
    stack->SetMaximum(minmax[1]);
  }
  else if( logy ){
    c->SetLogy();
    stack->SetMinimum(10e-2);
    stack->SetMaximum(minmax[1] * 100);
  }

  stack->Draw("hist");
  leg->Draw("same");
  TLatex* lat = makeLatex();
  lat->SetTextSize(0.06);
  for(uint l=0; l<lbls.size(); ++l)
    lat->DrawLatex( lbls[l].x, lbls[l].y, lbls[l].lbl.c_str() );
  
  c->SaveAs(save.c_str());
  
}

//---------------------------------------------------------------------//
// Build the THStack from a set of histograms
//---------------------------------------------------------------------//
THStack* FakePlotting::buildStack(TH1F* h[], int nhists)
{
  if(m_dbg) cout << "FakeClosurePlot::buildStack" << endl;

  int begin = nhists - 1; // Don't include total MC in the stack
  int end   = 0;                // Don't include data in the stack

  THStack* stack = new THStack("stack","stack");
  for(int i=begin; i>=end; i--){
    //cout<<"hist i: "<<i<<" "<<h[i]<<endl;
    stack->Add( h[i] );
  }
  return stack;

}

//--------------------------------------------------------//
// Get the Rate Histogram MC corrected
//--------------------------------------------------------//
TH1F* FakePlotting::getCorrectedRate(TFile* data, TFile* mc, string dataname, string mcname,
				     string xtitle, int color, int mark)
{
  cout<<"In Cor rate "<<data<<" "<<mc<<endl;
  TH1F* dtnum = Get(data, dataname+"_num", color, xtitle.c_str(), "Rate", mark);
  TH1F* dtden = Get(data, dataname+"_den", color, xtitle.c_str(), "Rate", mark);
  TH1F* mcnum = Get(mc,   mcname+"_num",   color, xtitle.c_str(), "Rate", mark);
  TH1F* mcden = Get(mc,   mcname+"_den",   color, xtitle.c_str(), "Rate", mark);
  mcnum->Scale(0.72);
  mcden->Scale(0.72);

  cout<<"Data den"<<endl;
  dumpHisto(dtden);
  cout<<"MC den"<<endl;
  dumpHisto(mcden);
  cout<<"Data num"<<endl;
  dumpHisto(dtnum);
  cout<<"MC num"<<endl;
  dumpHisto(mcnum);


  TH1F* correctedRate = (TH1F*) dtnum->Clone( (dataname+"_rat").c_str() );
  correctedRate->Reset();
  dtnum->Add(mcnum,-1.);
  dtden->Add(mcden,-1.);
  for(int bin=1; bin<=dtnum->GetNbinsX(); ++bin){
    if(dtnum->GetBinContent(bin) < 0) dtnum->SetBinContent(bin,0);
    if(dtden->GetBinContent(bin) < 0) dtden->SetBinContent(bin,0);
  }
  correctedRate->Divide(dtnum,dtden,1,1,"B");
  return correctedRate;
}

//--------------------------------------------------------//
TH1F* FakePlotting::buildLFRate(string lepton, string variable)				 
{

  string dataname = lepton + "_fakeLFZjet_" + variable;
  string mcname = lepton + "_fakeLFZjet_" + variable;

  TH1F* dtnum = Get(data.file, dataname+"_num", kBlack, "", "Rate", 20);
  TH1F* dtden = Get(data.file, dataname+"_den", kBlack, "", "Rate", 20);

  // Contamination
  TH1F* dibnum = Get(diboson.file, mcname+"_num", kBlack, "", "Rate", 20);
  TH1F* dibden = Get(diboson.file, mcname+"_den", kBlack, "", "Rate", 20);


  TH1F* correctedRate = (TH1F*) dtnum->Clone( (dataname+"_rat").c_str() );
  correctedRate->Reset();
  dtnum->Add(dibnum,-1.);
  dtden->Add(dibden,-1.);
  cout<<"dibnum: "<<dibnum->GetNbinsX()<<endl;
  cout<<"dibden: "<<dibden->GetNbinsX()<<endl;
  cout<<"datanum: "<<dtnum->GetNbinsX()<<endl;
  cout<<"dataden: "<<dtden->GetNbinsX()<<endl;

  for(int bin=1; bin<=dtnum->GetNbinsX(); ++bin){
    if(dtnum->GetBinContent(bin) < 0) dtnum->SetBinContent(bin,0);
    if(dtden->GetBinContent(bin) < 0) dtden->SetBinContent(bin,0);
  }
  correctedRate->Divide(dtnum,dtden,1,1,"B");
  return correctedRate;
}
//--------------------------------------------------------//
TH1F* FakePlotting::buildLFWjetRate(string lepton, string variable)
{

  string dataname = lepton + "_fakeLFWjet_" + variable;
  string mcname = lepton + "_fakeLFWjet_" + variable;

  TH1F* dtnum = Get(data.file, dataname+"_num", kBlack, "", "Rate", 20);
  TH1F* dtden = Get(data.file, dataname+"_den", kBlack, "", "Rate", 20);

  if(m_dbg) cout<<"dataname: "<<dataname<<"\nmcname: "<<mcname<<endl;

  // Contamination
  TH1F* dibnum   = Get(diboson.file, mcname+"_num", kBlack, "", "Rate", 20);
  TH1F* dibden   = Get(diboson.file, mcname+"_den", kBlack, "", "Rate", 20);
  cout<<"Have diboson"<<endl;
  TH1F* ttbarnum = Get(ttbar.file, mcname+"_num", kBlack, "", "Rate", 20);
  TH1F* ttbarden = Get(ttbar.file, mcname+"_den", kBlack, "", "Rate", 20);
  cout<<"Have ttbar"<<endl;
  TH1F* zjetnum  = Get(Zjets.file, mcname+"_num", kBlack, "", "Rate", 20);
  TH1F* zjetden  = Get(Zjets.file, mcname+"_den", kBlack, "", "Rate", 20);
  cout<<"Have zjets"<<endl;

  TH1F* correctedRate = (TH1F*) dtnum->Clone( (dataname+"_rat").c_str() );
  correctedRate->Reset();
  dtnum->Add(dibnum,-1.);
  dtden->Add(dibden,-1.);
  dtnum->Add(ttbarnum, -1);
  dtden->Add(ttbarden, -1);
  dtnum->Add(zjetnum, -1);
  dtden->Add(zjetden, -1);

  for(int bin=1; bin<=dtnum->GetNbinsX(); ++bin){
    if(dtnum->GetBinContent(bin) < 0) dtnum->SetBinContent(bin,0);
    if(dtden->GetBinContent(bin) < 0) dtden->SetBinContent(bin,0);
  }
  correctedRate->Divide(dtnum,dtden,1,1,"B");
  return correctedRate;
}

//--------------------------------------------------------//
// Build Rate histogram
//--------------------------------------------------------//
TH1F* FakePlotting::buildRate(TFile* file, string name, string sname, string xtitle, 
			      string ytitle, int color, int mark)
{
  TH1F* ratio=0;
  if(!file || !file->IsOpen()) {
    cout<<"buildRate: invalid input file"
        <<" ("<<file
        <<" : "<<(file ? file->GetName() : "")
        <<" "<<(file && file->IsOpen() ? "open":"not open")
        <<"), bailing out"<<endl;
    return ratio;
  }
  string hnameNum(name + "_num"), hnameDen(name + "_den");
  TH1F* num = Get(file, hnameNum, color, xtitle.c_str(), ytitle.c_str(), mark);
  TH1F* den = Get(file, hnameDen, color, xtitle.c_str(), ytitle.c_str(), mark);
  if(!num || !den) {
    cout<<"buildRate: invalid input histos"
        <<" "<<hnameNum<<"="<<num<<", "<<hnameDen<<"="<<den<<", bailing out"<<endl;
    return ratio;
  }  else if(m_dbg){
    cout<<"File ("<<file<<") : "<<file->GetName()<<endl;
    cout<<"Num: "<< hnameNum << endl;
    cout<<"Den: "<< hnameDen << endl;
  }

  //cout<<"\tnum: "<<num->GetMinimum()<<" den "<<den->GetMinimum()<<endl;
  //for(int bin=1; bin<=num->GetNbinsX(); ++bin)
  //cout<<"\t\tbin: "<<bin<<" "<<num->GetBinContent(bin)<<" "<<den->GetBinContent(bin)<<endl;
  ratio = (TH1F*) num->Clone( (name + "_rat").c_str() );
  ratio->Reset();
  ratio->Divide(num,den,1,1,"B");
  return ratio;
}

//--------------------------------------------------------//
// Build Side band subracted rate
//--------------------------------------------------------//
TH1F* FakePlotting::buildSideBandSubRate(TFile* file, string lep, string var, 
					 string xtitle, string ytitle, int color, int mark)
{

  // Sidebands:
  TH1F* numB = Get(file, lep+"_realSideLow_"+var+"_num",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* denB = Get(file, lep+"_realSideLow_"+var+"_den",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* numC = Get(file, lep+"_realSideHigh_"+var+"_num",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* denC = Get(file, lep+"_realSideHigh_"+var+"_den",color,xtitle.c_str(),ytitle.c_str(),mark);
  
  // Z window
  TH1F* numA = Get(file, lep+"_realCR_"+var+"_num",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* denA = Get(file, lep+"_realCR_"+var+"_den",color,xtitle.c_str(),ytitle.c_str(),mark);
  
  // Handle errors in a special way
  vector<float> errNum;
  vector<float> errDen;
  for(int bin = 1; bin<=numA->GetNbinsX(); ++bin){
    float errAnum = numA->GetBinError(bin);
    float errAden = denA->GetBinError(bin);
    float errBnum = numB->GetBinError(bin);
    float errBden = denB->GetBinError(bin);
    float errCnum = numC->GetBinError(bin);
    float errCden = denC->GetBinError(bin);
    errDen.push_back( sqrt(errAden*errAden - 
			   (errBden*errBden + errCden*errCden)) );
    errNum.push_back( sqrt(errAnum*errAnum - 
			   (errBnum*errBnum + errCnum*errCnum)) );

  }

  // Now have all histograms, and the errors for the bins.  Divide
  // and set the bin errors, then create ratio
  numB->Add(numC);
  denB->Add(denC);
  numA->Add(numB, -1);
  denA->Add(denB, -1);
  for(uint i=0; i<errNum.size(); ++i){
    numA->SetBinError(i+1, errNum.at(i));
    denA->SetBinError(i+1, errDen.at(i));
  }
  
  TH1F* ratio = (TH1F*) numA->Clone("rat");
  ratio->Reset();
  ratio->Divide(numA,denA,1,1,"B");

  numA->Delete();
  denA->Delete();  
  numB->Delete();
  denB->Delete();
  numC->Delete();
  denC->Delete();

  cout<<"FakePlotting::buildSideBandSubRate :"
      <<" ratio histogram for '"<<lep+"_"+var<<"'"<<endl;
  ratio->Print("all");
  return ratio;
  
}

//--------------------------------------------------------//
// Build Background Corrected Real
//--------------------------------------------------------//
TH1F* FakePlotting::buildCorrectedRealRate(TFile* file, string lep, string var, 
					   string xtitle, string ytitle, int color, int mark)
{

  // Need to take the rates from realCRB and subract them
  // from realCRA. CRB will require OS, so basically
  // taking the background values from em channel

  // Get background and signal
  TH1F* b_num = Get(file, lep+"_realCRD_actual_"+var+"_num",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* b_den = Get(file, lep+"_realCRD_actual_"+var+"_den",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* s_num = Get(file, lep+"_realCRA_actual_"+var+"_num",color,xtitle.c_str(),ytitle.c_str(),mark);
  TH1F* s_den = Get(file, lep+"_realCRA_actual_"+var+"_den",color,xtitle.c_str(),ytitle.c_str(),mark);

  //cout<<"Before: "<<b_num->Integral()<<endl;
  //b_num->Scale(0.5);
  //cout<<"After: "<<b_num->Integral()<<endl;
  //b_den->Scale(0.5);
  s_num->Add(b_num,-1);
  s_den->Add(b_den,-1);

  TH1F* ratio = (TH1F*) s_num->Clone("rat");
  ratio->Reset();
  ratio->Divide(s_num,s_den,1,1,"B");
  
  // Clean Up
  b_num->Delete();
  b_den->Delete();
  s_num->Delete();
  s_den->Delete();
  
  return ratio;

}

//--------------------------------------------------------//
// Build Relative Fraction for a region
//--------------------------------------------------------//
vector<float> FakePlotting::buildFraction(vector<File> files, string plotname)
{  

  string plot = plotname + "_den";
  vector<float> percentage;
  //percentage.resize( files.size() );
  float total = 0;

  for(uint f=0; f<files.size(); ++f){
    File F = files.at(f);
    if(F.sname == "mc") continue;
    percentage.push_back( ((TH1F*) F.file->Get(plot.c_str()))->Integral());
    total += percentage.back();
  }
  
  // Convert to percentage
  for(uint i=0; i<percentage.size(); ++i){
    if( total != 0 ) percentage[i] = (percentage[i]/total) * 100.;
    else percentage[i] = 0.;
  }


  return percentage;
}

//--------------------------------------------------------//
// Get Maximum
//--------------------------------------------------------//
float FakePlotting::getMaximum(TH1F* h[], uint nhisto)
{

  float max = 0;
  for(uint i=0; i<nhisto; ++i){
    int maxbin = h[i]->GetMaximumBin();
    float tempMax = h[i]->GetBinContent(maxbin) + h[i]->GetBinError(maxbin);
    if( max < tempMax )
      max = tempMax;
  }// end loop over histo
  
  return max;
}

//--------------------------------------------------------//
// Get Minimum
//--------------------------------------------------------//
float FakePlotting::getMinimum(TH1F* h[], uint nhisto)
{

  float min = 999;
  for(uint i=0; i<nhisto; ++i){
    int nbins = h[i]->GetNbinsX();
    for(int bin=1; bin<=nbins; ++bin){
      float tempMin = h[i]->GetBinContent(bin) - h[i]->GetBinError(bin);
      //if( tempMin < 0 ) tempMin = 0;
      if( min > tempMin )
	min = tempMin;
    }
  }// end loop over histo
  
  return min;
}

//--------------------------------------------------------//
// Dump histogram bin info
//--------------------------------------------------------//
void FakePlotting::dumpHisto(TH1F* h)
{
  h->Print();
  for(int bin = 1; bin<=h->GetNbinsX(); ++bin)
    cout<<"Bin: "<<bin<<" Content: "<<h->GetBinContent(bin)<<" +/- "<<h->GetBinError(bin)<<endl;

}

//--------------------------------------------------------//
// Check the percentages in each signal region
//--------------------------------------------------------//
void FakePlotting::checkPercentages()
{

  // Loop over the signal and control regions to try to 
  // determine what the statistical error on R_XR would
  // be.  If it is close for all SR, then we can assign
  // a conservative error to cover them all.  If it is
  // highly dependent, then will need to do one for each
  // or relax some cuts to let in more statistics.

  // Loop over lepton flavor
  for(int iL=0; iL<LT_N; ++iL){
    string lep = LTNames[iL];

    // What source to look at
    for(int iS=0; iS<LS_N; ++iS){
      if( iL == LT_MU && iS == LS_Conv ) continue;
      if( !(iS==LS_Conv || iS==LS_QCD || iS==LS_Real) ) continue;

      // Loop over control regions
      for(int iSR=0; iSR<SR_N; ++iSR){
	//if(iSR > VRSS) continue;
	string srname = SRNames[iSR];
    // DG 2013-09-16 this is really an ugly hack! to be removed ASAP.
    if(SR_WHSS==iSR) srname = CRNames[CR_SRWHSS];
	cout<<"--------------------------------------"<<endl;
	cout<<srname<<" for "<<LSNames[iS]<<endl;
      
	string pname = lep+"_"+srname+"_all_flavor_den";
	cout<<"\tPlot: "<<pname<<endl;

	vector<float> yields;
	vector<float> errs;
	float total = 0;
	// Loop over the files
	for(uint f = 0; f<m_files.size(); ++f){
	  File file = m_files.at(f);
	  TH1F* temp = (TH1F*) file.file->Get(pname.c_str());
      if(!temp) {
        cout<<"cannot get '"<<pname<<"' from '"<<file.file->GetName()<<"'...skip"<<endl;
        yields.push_back(0.0);
        errs.push_back(0.0);
        continue;
      }
	  yields.push_back(temp->GetBinContent(iS+1));
	  errs.push_back(temp->GetBinError(iS+1));
	  total += yields.back();
	  //cout<<"\t"<<file.name<<" = "<<temp->GetBinContent(iS+1)<<" +/- "<<temp->GetBinError(iS+1)<<endl;	 
	  
	}// end loop over files
	float totalErr = 0;
	for(uint f = 0; f<m_files.size(); ++f){
	  string name = m_files.at(f).name;
	  float yield = yields.at(f);
	  float weight = (total!=0.0 ? yield/total : 0.0);
	  cout<<"\t"<<name
          <<" Total: "<<yields.at(f)<<" +/- "<<errs.at(f)
          <<" weighted error: "<<weight*errs.at(f)<<endl;
	  float err = weight * (yield == 0 ? 0 : errs.at(f)/yield);
	  totalErr += err*err; 
	}// end loop over errs
	cout<<"\t\tTotal Error: "<<sqrt(totalErr)<<endl;
      }// end loop over SR
      
    }// end loop over Lepton Sources

  }// end loop over lepton types
}
//----------------------------------------------------------
bool FakePlotting::isCanvasWithNeededSf(const std::string &canvasName)
{
  vector<string> interestingCnames;
  interestingCnames.push_back("sf_elec_fakeConv_all_l_pt");
  interestingCnames.push_back("sf_elec_fakeHF_all_l_pt");
  interestingCnames.push_back("sf_elec_realCR_all_l_pt");
  interestingCnames.push_back("sf_muon_fakeHF_all_l_pt");
  interestingCnames.push_back("sf_muon_realCR_all_l_pt");
  vector<string>::const_iterator it=interestingCnames.begin(), end=interestingCnames.end();
  for(; it!=end; ++it) if(contains(canvasName, *it)) return true;
  return false;
}
//----------------------------------------------------------
