#ifndef fancyplotting_h
#define fancyplotting_h

//////////////////////////////////////////////
// Plotting script to obtain all formatted  //
// histograms for gamma+jet analysis.       //
//////////////////////////////////////////////

// Mine
#include "SusyTest0/myHist.h"
#include "SusyTest0/SusyPlotter.h"
#include "SusyTest0/SusyAnaDefsMatt.h"

// Root
#include "TF1.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

// Standard
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <assert.h>
#include <vector>
#include <sstream>

using namespace std;

// Enums
enum FPRunOption {
  RO_ALL = 0,
  RO_SR1,
  RO_SR2,
  RO_SR3,
  RO_SR4,
  RO_SR5,
  RO_ZWindow,
  RO_VR1,
  RO_VR2,
  RO_OSInc,
  RO_SSInc,
  RO_VR3,
  RO_VR4,
  RO_VRTL,
  RO_NONE,
  N_FPRunOption
};

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
  bool ismc;    // true if mc, false if data
  bool isfake;  // true if fake, false if not
  bool xsLumi;  // this is to be used if the stat error is zero
  bool isZ;     // to be used for Z+Met QCD
};

class FakeClosurePlot : public myHist {

 public:
  FakeClosurePlot();
  ~FakeClosurePlot();
  void init(FPRunOption opt = RO_ALL); // Initialize all the file objects
  void setPlots();
  string u(){ return "_"; };
  void DataMCAnaPlots();

  void buildHists(vector<TH1F*> &hists, vector<TH1F*> &sys, string var,
                  string xtitle, Chan ch, PlotRegion PR);
  TLegend* buildLegend(vector<TH1F*> hists, TGraphAsymmErrors* errs, float* x, float* y);
  TH1F* buildRatio(TH1F* data, TH1F* SM);
  THStack* buildStack(vector<TH1F*> hists);
  void plotAll(vector<TH1F*> hists, vector<TGraphAsymmErrors*> errs,
               string save, TLegend* leg, int ch, bool logy=false, bool logx=false);
  void clear(){  // Clear Hists
    for(uint i=0; i<m_hists.size(); ++i)
      if(m_hists.at(i)) m_hists.at(i)->Delete();
    m_hists.clear();
    for(uint i=0; i<m_sys.size(); ++i)
      if(m_sys.at(i)) m_sys.at(i)->Delete();
    m_sys.clear();
    for(uint i=0; i<m_errs.size(); ++i)
      if(m_errs.at(i)) m_errs.at(i)->Delete();
    m_errs.clear();

  }
  TGraphAsymmErrors* buildErrors(TH1F* summary, vector<TH1F*> sys);
  TGraphAsymmErrors* buildRatioErrors(TH1F* nominal, TGraphAsymmErrors* tg_errs);
  void addSysError(TH1F* nominal, TFile* file, string plot, vector<TH1F*> &sys);
  void addFakeSys(TH1F* nominal, TFile* file, string plot, vector<TH1F*> &sys);
  void getFakeSys(TH1F* nominal, TFile* file, string plot, float &sysup, float &sysdn);
  // Miscellaneous
  float getNorm(TH1* h);
  float getStat(TH1* h, float low, float high);
  void setMinMax(TH1* &h, float min, float max);
  void setMinMax(TH1F* &h, float min, float max){ setMinMax( (TH1*&) h, min, max ); };
  void setMinMax(TH1D* &h, float min, float max){ setMinMax( (TH1*&) h, min, max ); };
  float getMax(TH1F* h[], int n);
  TLine* getLine(TH1* h, float y, int color, int style);
  TLine* getLine(TH1F* h, float y, int color, int style){return getLine((TH1*)h,y,color,style);};
  TLine* getLine(TH1D* h, float y, int color, int style){return getLine((TH1*)h,y,color,style);};
  Label makeLabel(string l, float posx, float posy){Label l0; l0.lbl=l; l0.x=posx; l0.y=posy; return l0;};
  void setDebug(int d){ m_dbg = d; };
  void setIntegralOption(bool val){ m_addIntegral = val; };

  FakeClosurePlot& setTag(const std::string &name);
  FakeClosurePlot& setInputDir(const std::string &dir);
  FakeClosurePlot& setOuputDir(const std::string &dir);

 private:
  vector<File> m_files;                  // Vector for holding File objects
  vector< pair<string,string> > m_plots; // vector of plots
  vector<PlotRegion> m_PRs;              // Plot regions to plot
  vector<TH1F*> m_hists;                 // vector of histograms
  vector<TH1F*> m_sys;                   // vector for sys histograms -- 2 entries
  vector<TGraphAsymmErrors*> m_errs;     // Error objects
  FPRunOption m_opt;                     // Store run option
  int m_dbg;                             // Debug flag
  bool m_addIntegral;                    // Add integral to legend
  int m_MCColor;                         // Color for the total MC -- kRed
  string m_tag;                          // tag of the input files
  string m_inputdir;                     // dir of the input files
  string m_outputdir;
};

#endif
