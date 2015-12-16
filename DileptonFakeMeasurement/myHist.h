#ifndef MYHIST_H
#define MYHIST_H


//-------------------------------------------------------//
// This class is a nice way to add more methods to 
// manipulate histograms in a clean way for the analysis
// class.  At least that is it's initial purpose...
//-------------------------------------------------------//

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TLine.h"
#include "THStack.h"
#include "TF1.h"

#include <string>
#include <iostream>


using namespace std;

class myHist{

 public:
  myHist(){ SetStyle(); };
  virtual ~myHist(){};
  
  //------------------------//
  // Create 1-D Histograms
  //------------------------//
  TH1F* hist(string name, string xtitle, const int nbins, const float* bins){    
    TH1F* h = new TH1F(name.c_str(),name.c_str(),nbins,bins);
    h->GetXaxis()->SetTitle( xtitle.c_str() );
    h->Sumw2();
    return h;
  };
  
  TH1F* hist(string name, string xtitle, const int nbins, const float min, const float max){    
    TH1F* h = new TH1F(name.c_str(),name.c_str(),nbins,min,max);
    h->GetXaxis()->SetTitle( xtitle.c_str() );
    h->Sumw2();
    return h;
  };

  //------------------------//
  // Create 2-D Histograms
  //------------------------//
  TH2F* hist2(string name, string xtitle, string ytitle, 
	      const int nxbins, const float* xbins,
	      const int nybins, const float* ybins)
    {    
      TH2F* h = new TH2F(name.c_str(),name.c_str(),nxbins,xbins,nybins,ybins);
      h->GetXaxis()->SetTitle( xtitle.c_str() );
      h->GetYaxis()->SetTitle( ytitle.c_str() );
      h->Sumw2();
      return h;
    };
  TH2F* hist2(string name, string xtitle, string ytitle, 
	      const int nxbins, const float* xbins,
	      const int nybins, const float ymin, const float ymax)
    {    
      float ybins[nybins];
      float delta = (ymax-ymin)/nybins;
      for(int i=0; i<=nybins; ++i)
	ybins[i] = ymin + i*delta;

      return hist2(name,xtitle,ytitle,nxbins,xbins,nybins,ybins);
    };
  TH2F* hist2(string name, string xtitle, string ytitle, 
	      const int nxbins, const float xmin, const float xmax,
	      const int nybins, const float* ybins)
    {    
      float xbins[nxbins];
      float delta = (xmax-xmin)/nxbins;
      for(int i=0; i<=nxbins; ++i)
	xbins[i] = xmin + i*delta;

      return hist2(name,xtitle,ytitle,nxbins,xbins,nybins,ybins);
    };
  TH2F* hist2(string name, string xtitle, string ytitle, 
	      const int nxbins, const float xmin, const float xmax,
	      const int nybins, const float ymin, const float ymax)
    {
      TH2F* h = new TH2F(name.c_str(),name.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
      h->GetXaxis()->SetTitle( xtitle.c_str() );
      h->GetYaxis()->SetTitle( ytitle.c_str() );
      h->Sumw2();
      return h;
    };

  //-----------------------------------//
  // Make TProfile
  //-----------------------------------//
  TProfile* prof(string name, const int nxbins, const float* xbins){

    TProfile* pro = new TProfile(name.c_str(), name.c_str(), nxbins, xbins);
    pro->Sumw2();
    return pro;

  };

  TProfile* prof(string name, const int nxbins, const float xmin, const float xmax){

    TProfile* pro = new TProfile(name.c_str(), name.c_str(), nxbins, xmin, xmax);
    pro->Sumw2();
    return pro;

  };

  //-----------------------------------//
  // Format histogram
  //-----------------------------------//
  void setAttributes(TH1F* &h, int color, int marker=20){
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
  };
  void setAttributes(TH1D* &h, int color, int marker=20){
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
  };

  //-----------------------------------//
  // Make useful plot objects
  //-----------------------------------//
  TLine* makeLine(float x0, float x1, float y0, float y1, int color=kBlack, int style=1){
    TLine* line = new TLine(x0,y0,x1,y1);
    line->SetLineWidth(1);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    return line;
  };

  TLegend* makeLegend(float* x, float* y){

    TLegend* l = new TLegend(x[0],y[0],x[1],y[1]);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetLineColor(0);
    l->SetBorderSize(0);
    l->SetTextSize( 0.05 );

    return l;
  };

  TLatex* makeLatex(){
    TLatex* lat = new TLatex();
    lat->SetTextSize(0.07);
    lat->SetNDC();
    lat->SetTextFont(42);
    lat->SetTextColor(kBlack);
    return lat;
  };

  TCanvas* makeCanvas(string name){
    TCanvas* c_new = new TCanvas(name.c_str(), name.c_str(), 600, 500);
    c_new->SetTopMargin(0.05);
    c_new->SetBottomMargin(0.12);
    c_new->SetLeftMargin(0.12);
    c_new->SetRightMargin(0.05);
    return c_new;
  };

  //-----------------------------------//
  // Lumi and sqrt 
  //-----------------------------------//
  void AddLumi(TVirtualPad* &tv, string lumi, float x, float y){

    tv->cd();
    
    TLatex l;
    l.SetTextSize(0.04);
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextColor(kBlack);
    l.DrawLatex(x,y,Form("#int L dt = %s",lumi.c_str()));
    l.DrawLatex(x,y-0.1,"#sqrt{s} = 7 TeV");

  };

  //-----------------------------------//
  // Functions to retrieve Histograms
  //-----------------------------------//
  TH1F* Get(TFile* file, string name, int color=kBlack, const char* title ="",const char* ytitle="",int marker=20){
    TH1F *h=0;
    if(!file) {
      cout<<"Get: invalid input file "<<file<<endl;
      return h;
    }
    h = (TH1F*) file->Get(name.c_str());
    if(!h) {
      cout<<"cannot get histo '"<<name<<"' from '"<<file->GetName()<<"'"<<endl;
      return h;
    }
    if(sizeof(title) != 0)  h->GetXaxis()->SetTitle(title);
    h->GetYaxis()->SetTitle(ytitle);
    //h->SetTitle(title);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.75);
    //h->Scale(1/h->Integral());
    return h;
  };
  
  TH1D* GetP(TFile* file, string name, string projname, int color=kBlack, 
	     const char* xtitle = "", const char* ytitle = "", int marker=20){
    TProfile* p = (TProfile*) file->Get(name.c_str());
    TH1D* h = (TH1D*) p->ProjectionX(projname.c_str());
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    setAttributes(h, color, marker);
    return h;
  };

  TProfile* GetProf(TFile* file, string name){
    return (TProfile*) file->Get(name.c_str());
  };

  TH2F* Get2(TFile* file, string name){
    return (TH2F*) file->Get(name.c_str());
  };

  //------------------------------------//
  // Make a ratio histogram
  //------------------------------------//
  TH1F* RatioHist(TH1* h_top, TH1* h_bot, const char* label, int color=kBlack, string opt = ""){

    // Make the ratio plot:
    TH1F* ratio = (TH1F*) h_top->Clone("ratio");
    ratio->Reset();
    ratio->GetXaxis()->SetTitle(h_top->GetXaxis()->GetTitle());
    ratio->GetYaxis()->SetTitle(label);
    ratio->Divide(h_top,h_bot,1,1,opt.c_str());
    ratio->SetLabelSize(0.12,"X");
    ratio->SetLabelSize(0.12,"Y");
    ratio->SetTitleFont(42,"X");
    ratio->SetTitleFont(42,"Y");
    ratio->GetXaxis()->SetTitleSize(0.15);
    ratio->GetXaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetNdivisions(205);
    ratio->GetYaxis()->SetTitleSize(0.11);
    ratio->GetYaxis()->SetTitleOffset(0.55);
    ratio->SetLineColor(color);
    ratio->SetMarkerColor(color);

    return ratio;
  };

  //------------------------------------//
  // Create TVirtual Pad for ratio plots
  //------------------------------------//
  TVirtualPad* CreatePad(TCanvas* c, TPad* &_pTop, TPad* &_pBot, float ysep = 0.3){
    
    // Create Virtual Pad
    TVirtualPad* _tv = c->cd();
    _pTop = new TPad("pTop","pTop",0,ysep,1,1);
    _pTop->SetTopMargin(0.05);
    _pTop->SetBottomMargin(0.0);
    _pTop->SetRightMargin(0.05);
    _pTop->SetLeftMargin(0.12);
    _pTop->SetNumber(1);

    _pBot = new TPad("pBot","pBot",0,0,1,ysep);
    _pBot->SetTopMargin(0.00);
    _pBot->SetBottomMargin(0.25);
    _pBot->SetRightMargin(0.05);
    _pBot->SetLeftMargin(0.12);
    _pBot->SetNumber(2);
    _tv->cd();

    return _tv;
  };

  //------------------------------------//
  // Moving away from the virtual pad
  //------------------------------------//
  void makePads(TCanvas* &c, TPad* &_pTop, TPad* &_pBot, float ysep = 0.3){
    
    c->cd();
    _pTop = new TPad("pTop","pTop",0,ysep,1,1);
    _pTop->SetTopMargin(0.05);
    _pTop->SetBottomMargin(0.0);
    _pTop->SetRightMargin(0.05);
    _pTop->SetLeftMargin(0.12);
    _pTop->SetNumber(1);

    _pBot = new TPad("pBot","pBot",0,0,1,ysep);
    _pBot->SetTopMargin(0.00);
    _pBot->SetBottomMargin(0.3);
    _pBot->SetRightMargin(0.05);
    _pBot->SetLeftMargin(0.12);
    _pBot->SetNumber(2);

  };

  //------------------------------------//
  // Get max of two histograms
  //------------------------------------//
  double getMax(TH1* h0, TH1* h1){

    double max0=(h0->GetMaximum()+h0->GetBinError(h0->GetMaximumBin()));
    double max1=(h1->GetMaximum()+h1->GetBinError(h1->GetMaximumBin()));
    return (max0<max1) ? max1 : max0;

  };

  //------------------------------------//
  // Set Min max based on histograms
  //------------------------------------//
  vector<float> getMinMax(TH1F* h[], unsigned int nhist){

    vector<float> minmax;

    if(nhist == 0) return minmax;
    float min = 999;
    float max = -1; 
    for(unsigned int i=0; i<nhist; ++i){
      for(int bin = 1; bin <= h[i]->GetNbinsX(); ++bin){
	float bc = h[i]->GetBinContent(bin);
	float be = h[i]->GetBinError(bin);
	//cout<<"Min: "<<min<<" max: "<<max<<" bc: "<<bc<<" +/- "<<be<<endl;
	if(bc == 0) continue;
	if( min > bc - be )
	  min = bc - be;
	if( max < bc + be )
	  max = bc + be;

      }// end loop over bins
    }// end loop over hists

    if(min < 0) min = 0;
    minmax.push_back(min);
    minmax.push_back(max);
    return minmax;
  };

  // from a vector of hists
  vector<float> getMinMax(vector<TH1F*> h){

    vector<float> minmax;
    unsigned int nhist = h.size();

    if(nhist == 0) return minmax;
    float min = 999;
    float max = -1; 
    for(unsigned int i=0; i<nhist; ++i){
      for(int bin = 1; bin <= h[i]->GetNbinsX(); ++bin){
	float bc = h[i]->GetBinContent(bin);
	float be = h[i]->GetBinError(bin);
	//cout<<"Min: "<<min<<" max: "<<max<<" bc: "<<bc<<" +/- "<<be<<endl;
	if(bc == 0) continue;
	if( min > bc - be )
	  min = bc - be;
	if( max < bc + be )
	  max = bc + be;

      }// end loop over bins
    }// end loop over hists

    if(min < 0) min = 0;
    minmax.push_back(min);
    minmax.push_back(max);
    return minmax;
  };

  //-----------------------------------//
  // Set Custom Style
  //-----------------------------------//
  TStyle* atlasStyle;
  void SetStyle(){

    atlasStyle = new TStyle("ATLAS","Atlas style");
    Int_t icol=0;
    Int_t font = 42;
    Double_t tsize=0.05;
    Double_t tlabelsize=0.04;
    Int_t markerStyle=20;
    Double_t msize =1.0; 
    Int_t hlinewidth =1.0;

    // Canvas settings   
    atlasStyle->SetFrameBorderMode(icol);
    atlasStyle->SetFrameFillColor(icol);
    atlasStyle->SetCanvasBorderMode(icol);
    atlasStyle->SetPadBorderMode(icol);
    atlasStyle->SetPadColor(icol);
    atlasStyle->SetCanvasColor(icol);
    atlasStyle->SetStatColor(icol);

    // set the paper & margin sizes
    atlasStyle->SetPaperSize(20,26);
    atlasStyle->SetPadTopMargin(0.05);
    atlasStyle->SetPadRightMargin(0.12);
    atlasStyle->SetPadBottomMargin(0.15);
    atlasStyle->SetPadLeftMargin(0.18);
    
    // use large fonts                 
    atlasStyle->SetTextFont(font);
    atlasStyle->SetTextSize(tsize);

    atlasStyle->SetTitleFont(font,"x");
    atlasStyle->SetTitleFont(font,"y");
    atlasStyle->SetTitleFont(font,"z");
    atlasStyle->SetTitleOffset(1.4,"x");
    atlasStyle->SetTitleOffset(1.5,"y");
    atlasStyle->SetTitleOffset(1.,"z");
    atlasStyle->SetTitleSize(tsize,"x");
    atlasStyle->SetTitleSize(tsize,"y");
    atlasStyle->SetTitleSize(tsize,"z");
    atlasStyle->SetTickLength (0.02,"x");
    atlasStyle->SetTickLength (0.02,"y");
    atlasStyle->SetTickLength (0.02,"z");

    atlasStyle->SetLabelFont(font,"x");
    atlasStyle->SetLabelFont(font,"y");
    atlasStyle->SetLabelFont(font,"z");
    //  atlasStyle->SetLabelOffset( 0.01, "x");
    atlasStyle->SetLabelOffset( 0.02, "y");
    //  atlasStyle->SetLabelOffset( 0.01, "z");
    atlasStyle->SetLabelSize(tlabelsize,"x");
    atlasStyle->SetLabelSize(tlabelsize,"y");
    atlasStyle->SetLabelSize(tlabelsize,"z");

    //palette settings                         
    atlasStyle->SetPalette(1,0);

    //use bold lines and markers               
    atlasStyle->SetMarkerStyle(markerStyle);
    atlasStyle->SetMarkerSize(msize);
    atlasStyle->SetHistLineWidth(hlinewidth);
    atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    atlasStyle->SetEndErrorSize(0.); //remove erro bar caps          

    //do not display any of the standard histogram decorations       
    atlasStyle->SetStatX(0.99);
    atlasStyle->SetStatY(0.99);
    atlasStyle->SetStatH(0.01);
    atlasStyle->SetStatW(0.2);
    
    atlasStyle->SetStatStyle(0);
    atlasStyle->SetStatFont(font);
    atlasStyle->SetStatFontSize(0.03);
    atlasStyle->SetOptStat(0);
    atlasStyle->SetStatBorderSize(1);
    atlasStyle->SetOptTitle(0);
    atlasStyle->SetOptFit(0);

    atlasStyle->SetTitleStyle(icol);
    atlasStyle->SetTitleH(0.08);

    // put tick marks on top and RHS of plots
    atlasStyle->SetPadTickX(1);
    atlasStyle->SetPadTickY(1);

    gROOT->SetStyle("ATLAS");

  };

  ClassDef(myHist, 1);

};

#endif
