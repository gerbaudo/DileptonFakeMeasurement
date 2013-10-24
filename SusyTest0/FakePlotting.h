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
  RO_DataMCSF,
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
  void DataMCSF();

  float getBC(TFile* f, string plot, int bin=1){
    return ((TH1F*) f->Get(plot.c_str()))->GetBinContent(bin);
  };
  float getBE(TFile* f, string plot, int bin=1){
    return ((TH1F*) f->Get(plot.c_str()))->GetBinError(bin);
  };
  void plotSF(TH1F* h[], vector<string> names, TCanvas* c,
              float* MinMax, float* xLeg, float* yLeg,
              vector<Label> lbls, string save, bool doFit = true,
              bool topBot = true /*data/MC = true, MC/Data=false*/,
              string topLabel = "Rate", string ratLabel = "Data/MC");
  TH1F* buildRate(TFile* file, string name, string sname, string xtitle,
                  string ytitle, int color, int mark=20);

  // Build RE histogram doing sideband subtraction
  TH1F* buildSideBandSubRate(TFile* file, string lep, string var, string xtitle,
                             string ytitle, int color, int mark=20);
  vector<float> buildFraction(vector<File> files, string plotname);
  float getMaximum(TH1F* h[], uint nhisto);
  float getMinimum(TH1F* h[], uint nhisto);
  void setDebug(int dbg){ m_dbg = dbg; };
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
