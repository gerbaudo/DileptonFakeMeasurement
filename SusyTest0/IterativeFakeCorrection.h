#ifndef IterativeFakeCorrection_h
#define IterativeFakeCorrection_h

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "SusyTest0/myHist.h"
#include "SusyTest0/SusyAnaDefsMatt.h"

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

class IterativeFakeCorrection : public myHist
{

 public:
  IterativeFakeCorrection();
  virtual ~IterativeFakeCorrection();

  //
  // Methods to perform iterations
  //

  void iterate();
  void correctRate(TH1F* &rate, TH1F* data, TH1F* mc, vector<double> corrections);
  vector<double> getC(TH1F* rate, TH1F* data_num, TH1F* data_den,  TH1F* mc, bool tight);

  
  // 
  // Methods to retrieve and modify histograms
  //

  TH1F* getHist(TFile* file, string histname);

  // 
  // Methods to get the weights
  //
  
  double getFake(float r, float f, float nLoose, float nTight, bool tight);

  void setDebug(bool dbg){ m_dbg = dbg; };
  void setNIter(int niter){ m_nIter = niter; };
  IterativeFakeCorrection& setInputMc(const std::string &filename);
  IterativeFakeCorrection& setInputData(const std::string &filename);
  IterativeFakeCorrection& setOutputFilename(const std::string &filename);
 protected:
  
  bool m_dbg;

  TH1F* m_real;
  File m_data;
  File m_mc;
  int m_nIter;
  std::string m_outputFilename;
};

#endif
