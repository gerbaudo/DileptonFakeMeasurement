#ifndef EffObject_h
#define EffObject_h

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// For some reason TEfficiency Objects are not able   //
// to be written to TFile's and I am tired of trying  //
// to use them.  So I will stick with histograms, but //
// make use of some convenience using this new        //
// defined data structure : )                         //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "TH1F.h"
#include "TH2F.h"
//#include "TEfficiency.h"

class EffObject : public TObject
{
 public:

  //-----------------------------------------//
  // Constructors
  //-----------------------------------------//
  EffObject(string name, const int nbins, const float* bins){
    num = new TH1F((name+"_num").c_str(),name.c_str(),nbins,bins);
    den = new TH1F((name+"_den").c_str(),name.c_str(),nbins,bins);
    num->Sumw2();
    den->Sumw2();
  };
  EffObject(string name, const int nbins, const float min, const float max){
    num = new TH1F((name+"_num").c_str(),name.c_str(),nbins,min,max);
    den = new TH1F((name+"_den").c_str(),name.c_str(),nbins,min,max);
    num->Sumw2();
    den->Sumw2();
  };

  //-----------------------------------------//
  // Destructor
  //-----------------------------------------//
  ~EffObject(){
    if(num) num->Delete();
    if(den) den->Delete();
  };

  void Fill(bool pass, float weight, float var){
    
    // Don't go outside bounds
    float max = den->GetXaxis()->GetXmax();			
    float min = den->GetXaxis()->GetXmin();
    float fillx = var > max ? max - 1e-4 : var;	
    fillx = fillx < min ? min + 1e-4 : fillx;

    den->Fill(fillx,weight);
    if(pass) num->Fill(fillx,weight);
  };

  void SetXLabel(int bin, string label){
    this->num->GetXaxis()->SetBinLabel(bin, label.c_str());
    this->den->GetXaxis()->SetBinLabel(bin, label.c_str());
  }

  ClassDef(EffObject,1);

 private:
  TH1F* num;
  TH1F* den;
  
};

// 2D Eff object
class EffObject2 : public TObject
{
 public:

  //-----------------------------------------//
  // Constructors
  //-----------------------------------------//
  EffObject2(string name, const int nxbins, const float* xbins,
	     const int nybins, const float* ybins){
    num = new TH2F((name+"_num").c_str(),name.c_str(),nxbins,xbins,nybins,ybins);
    den = new TH2F((name+"_den").c_str(),name.c_str(),nxbins,xbins,nybins,ybins);
    num->Sumw2();
    den->Sumw2();
  };
  EffObject2(string name, const int nxbins, const float xmin, const float xmax,
	     const int nybins, const float ymin, const float ymax){
    num = new TH2F((name+"_num").c_str(),name.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
    den = new TH2F((name+"_den").c_str(),name.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
    num->Sumw2();
    den->Sumw2();
  };

  //-----------------------------------------//
  // Destructor
  //-----------------------------------------//
  ~EffObject2(){
    if(num) num->Delete();
    if(den) den->Delete();
  };

  void Fill(bool pass, float weight, float varx, float vary){
    
    // Don't go outside bounds
    float xmax = den->GetXaxis()->GetXmax();			
    float xmin = den->GetXaxis()->GetXmin();
    float ymax = den->GetYaxis()->GetXmax();			
    float ymin = den->GetYaxis()->GetXmin();
    float fillx = varx > xmax ? xmax - 1e-4 : varx;	
    float filly = vary > ymax ? ymax - 1e-4 : vary;	
    fillx = fillx < xmin ? xmin + 1e-4 : fillx;
    filly = filly < ymin ? ymin + 1e-4 : filly;

    den->Fill(fillx,filly,weight);
    if(pass) num->Fill(fillx,filly,weight);
  };

  void SetXLabel(int bin, string label){
    this->num->GetXaxis()->SetBinLabel(bin, label.c_str());
    this->den->GetXaxis()->SetBinLabel(bin, label.c_str());
  }
  void SetYLabel(int bin, string label){
    this->num->GetYaxis()->SetBinLabel(bin, label.c_str());
    this->den->GetYaxis()->SetBinLabel(bin, label.c_str());
  }

  ClassDef(EffObject2,1);

 private:
  TH2F* num;
  TH2F* den;
  
};

#endif
