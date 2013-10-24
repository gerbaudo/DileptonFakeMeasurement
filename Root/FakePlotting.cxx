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
  m_runopt(runOpt),
  m_fakeStat(NULL)
{
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
  if(!inputIsValid) cout<<"invalid input '"<<fname<<"'"<<" ("<<out.file<<")"<<endl;
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
  m_outputdir = mkdirIfNeeded(endswith(dir,"/") ? dir : dir+"/");
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
  string tag = "_"+m_tag;
  initInputFile(data,    inDir+"data"       +tag+".root", "Data",        "data",  kBlack,    20);
  initInputFile(totMC,   inDir+"allBkg"     +tag+".root", "MC Combined", "mc",    kRed,      25);
  initInputFile(ttbar,   inDir+"ttbar"      +tag+".root", "t#bar{t}",    "ttbar", kRed,      21);
  initInputFile(Zjets,   inDir+"zjets"      +tag+".root", "Z+jets",      "Zjets", kOrange,   23);
  initInputFile(Wjets,   inDir+"wjets"      +tag+".root", "W+jets",      "Wjets", kViolet+2, 22);
  initInputFile(diboson, inDir+"diboson"    +tag+".root", "Diboson",     "dib",   kMagenta,  34);
  initInputFile(HF,      inDir+"heavyflavor"+tag+".root", "HF",          "HF",    kBlue,     26);

  if(m_runopt == RO_DataMCSF){
    m_files.push_back(data);
    m_files.push_back(totMC);
  }
}

//--------------------------------------------------------//
// Destructor
//--------------------------------------------------------//
FakePlotting::~FakePlotting()
{

  // Clean up the files
  for(uint i=0; i<m_files.size(); ++i) m_files.at(i).file->Delete();
  m_files.clear();
  delete m_fakeStat;
}
//--------------------------------------------------------//
// Get SF for Data and MC
//--------------------------------------------------------//
void FakePlotting::DataMCSF()
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
// Get Maximum
//--------------------------------------------------------//
float FakePlotting::getMaximum(TH1F* h[], uint nhisto)
{
  float max = 0;
  for(uint i=0; i<nhisto; ++i){
    int maxbin = h[i]->GetMaximumBin();
    float tempMax = h[i]->GetBinContent(maxbin) + h[i]->GetBinError(maxbin);
    if( max < tempMax ) max = tempMax;
  }
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
      if( min > tempMin ) min = tempMin;
    }
  }
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
