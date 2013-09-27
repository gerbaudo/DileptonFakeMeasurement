#ifndef FakePlotting_h
#define FakePlotting_h


////////////////////////////////////////////////////
// Script for formatting the fake rate and real   //
// efficincy plots.  All control regions in data  //
// and MC are handled.  Also this code will spit  //
// out the relative contributions from various    //
// fake 'flavors' in the differnt control regions //
////////////////////////////////////////////////////


#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "SusyTest0/myHist.h"
#include "SusyTest0/SusyAnaDefsMatt.h"
#include "SusyMatrixMethod/StatErrorTool.h"

#include "TMath.h"


using namespace std;

// Some generic structs to help with plotting
struct Label {
  string lbl;
  float x;
  float y;
};

struct File {
  TFile* file;  // Pointer to file
  string name;  // Name for plots, "Data 2012" "Pythia"
  string sname; // short name, like py, or dt11
  int color;    // color for the histogram
  int marker;   // maker style
};

// Typedefs
typedef unsigned int uint;

enum RunOption
{
  RO_Data = 0,
  RO_MC,
  RO_DataMC,
  RO_DataMCSub,
  RO_GJetCR,
  RO_DataMCSF,
  RO_SRRates,
  RO_SRComp,
  RO_TLInfo,
  RO_SRDump,
  RO_HFNorm,
  RO_DumpPer,
  RO_N
};

class FakePlotting : public myHist
{

 public:

  // Contstructor/Destructor
  FakePlotting(RunOption);
  virtual ~FakePlotting();
  void init();
  // initialize the File attributes; return false if input file is not there
  static bool initInputFile(File &out,
                            const string &fname, const string &name, const string &sname,
                            int color, int marker);
  FakePlotting& setTag(const std::string &name);
  FakePlotting& setInputDir(const std::string &dir);
  FakePlotting& setOuputDir(const std::string &dir);
  FakePlotting& setInputItercorrFile(const std::string &name);
  // The various looping functions
  void MCFakeRate();
  void DataFakeRate();
  void DataRateMCSub();
  void GammaJetCRRates();
  void DataMCSF(RunOption ro);
  void MCSRRates();
  void Composition();
  void TLPlots();
  void dumpSRTable();
  void GetHFNorm();

  // Methods to dump signal region fake
  void dumpSRFake();
  float getBC(TFile* f, string plot, int bin=1){
    return ((TH1F*) f->Get(plot.c_str()))->GetBinContent(bin);
  };
  float getBE(TFile* f, string plot, int bin=1){
    return ((TH1F*) f->Get(plot.c_str()))->GetBinError(bin);
  };
  float getNewStat(TFile* f, string sr, string ch);

  // Checking for the stat error
  void checkPercentages();

  // Plotting functions
  void plotFake(TH1F* h[], vector<string> names, TCanvas* c, 
		float* MinMax, float* xLeg, float* yLeg,
		vector<Label> lbls, string save, bool restrictMax=true,
		bool logy = false);
  void plotSF(TH1F* h[], vector<string> names, TCanvas* c, 
	      float* MinMax, float* xLeg, float* yLeg,
	      vector<Label> lbls, string save, bool doFit = true,
	      bool topBot = true /*data/MC = true, MC/Data=false*/,
	      string topLabel = "Rate", string ratLabel = "Data/MC");

  void plotStack(TH1F* h[], vector<string> names, TCanvas* c, 
		 float* xLeg, float* yLeg,
		 vector<Label> lbls, string save, bool logy = false);

  // Build Fake histogram
  TH1F* buildRate(TFile* file, string name, string sname, string xtitle, 
		  string ytitle, int color, int mark=20);
  TH1F* buildLFRate(string lepton, string variable);
  TH1F* buildLFWjetRate(string lepton, string variable);

  // Build RE histogram doing sideband subtraction
  TH1F* buildSideBandSubRate(TFile* file, string lep, string var, string xtitle, 
			     string ytitle, int color, int mark=20);
  TH1F* buildCorrectedRealRate(TFile* file, string lep, string var, string xtitle, 
			       string ytitle, int color, int mark=20);

  // Build Fractions
  vector<float> buildFraction(vector<File> files, string plotname);

  // Get Corrected Fake Rate
  TH1F* getCorrectedRate(TFile* data, TFile* mc, string dataname, string mcname,
			 string xtitle, int color, int mark = 20);

  // Wrapper to make a lebel struct
  Label makeLabel(string name, float xpos, float ypos){
    Label l;
    l.lbl = name;
    l.x = xpos;
    l.y = ypos;
    return l;
  };

  // Make stacked histogram
  THStack* buildStack(TH1F* h[], int nhists);

  // Set the maximum for histogram
  float getMaximum(TH1F* h[], uint nhisto);
  float getMinimum(TH1F* h[], uint nhisto);
  
  // Flag for debuggin'
  void setDebug(int dbg){ m_dbg = dbg; };

  // Dump contents of a histogram
  void dumpHisto(TH1F* h);
  static bool isCanvasWithNeededSf(const std::string &canvasName);
 protected:
  
  vector<File> m_files;

  File data;
  File totMC;
  File ttbar;
  File Zjets;
  File Ztandp;
  File Wjets;
  File stop;
  File diboson;
  File gjet;
  File fakePred;
  File HF;

  int m_dbg;
  
  RunOption m_runopt;

  FakeStatTool::StatErrorTool* m_fakeStat;
  std::string m_tag;
  std::string m_inputdir;
  std::string m_outputdir;
  std::string m_inputItercorrFile;

};









#endif
