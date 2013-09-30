// This plotting script will grow and be used to format
// histograms that can be shown in talks

#include "SusyTest0/FakeClosurePlot.h"
#include "SusyTest0/utils.h"

//---------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------//
FakeClosurePlot::FakeClosurePlot(/*FPRunOption opt*/) :
  myHist(),
  m_opt(RO_ALL),
  m_dbg(0),
  m_addIntegral(false),
  m_MCColor(kRed)
{
}

//---------------------------------------------------------------------//
// Initialize all files and things
//---------------------------------------------------------------------//
void initInputFile(File &out,
                   const string &fname, const string &name, const string &sname,
                   int color, bool ismc, bool isfake)
{
  out.file  = new TFile(fname.c_str());
  out.name  = name;
  out.sname = sname;
  out.color = color;
  out.ismc  = ismc;
  out.isfake = isfake;
  bool inputIsValid(out.file && out.file->IsOpen());
  if(!inputIsValid) cout<<"invalid input '"<<fname<<"'"<<" ("<<out.file<<")"<<endl;
}
//----------------------------------------------------------
void FakeClosurePlot::init(FPRunOption opt)
{
  m_opt = opt;
  File data, top, Zjets, Wjets, dib, hf, qcd;
  const string inDir(m_inputdir), tag(m_tag);
  const string tagg(tag+".AnaHists.root"), tagf(tag+".FakeHists.root");
  initInputFile(data,  inDir+"data"       +tagg, "Data",    "data",    kBlack,     false, false);
  initInputFile(top,   inDir+"ttbar"      +tagg, "Ttbar",   "ttbar",   kBlue,      true,  false);
  initInputFile(Zjets, inDir+"zjets"      +tagg, "Z+jet",   "Zjet",    kRed,       true,  false);
  initInputFile(Wjets, inDir+"wjets"      +tagg, "W+jet",   "Wjet",    kMagenta,   true,  false);
  initInputFile(dib,   inDir+"diboson"    +tagg, "Diboson", "diboson", kOrange,    true,  false);
  initInputFile(hf,    inDir+"heavyflavor"+tagg, "b-bbar",  "bbbar",   kSpring+1,  true,  false);
  initInputFile(qcd,   inDir+"data"       +tagf, "Fake",    "fake",    kGray,      true,  true );
  m_files.clear();

  m_files.push_back(data);
  m_files.push_back(top);
  m_files.push_back(Zjets);
  m_files.push_back(Wjets);
  m_files.push_back(dib);
  m_files.push_back(hf);
  m_files.push_back(qcd);

  if(opt == RO_ALL){ m_PRs.push_back(PR_VRSS); }
  else return;
  setPlots();
}
//---------------------------------------------------------------------//
// Destructor
//---------------------------------------------------------------------//
FakeClosurePlot::~FakeClosurePlot()
{

}
//---------------------------------------------------------------------//
// Initialize vectors
//---------------------------------------------------------------------//
void FakeClosurePlot::setPlots()
{
  // Set any and all plots that are needed and the functions will loop over this list
  m_plots.push_back( pair<string,string> ("l0_pt", "l_{0} P_{T} [GeV]") );
  m_plots.push_back( pair<string,string> ("l1_pt", "l_{1} P_{T} [GeV]") );
  m_plots.push_back( pair<string,string> ("ll_M", "m(ll) [GeV]") );
  m_plots.push_back( pair<string,string> ("met", "#slash{E}_{T} [GeV]") );
  m_plots.push_back( pair<string,string> ("metrel", "#slash{E}^{rel}_{T} [GeV]") );
  m_plots.push_back( pair<string,string> ("njets", "# jets") );
  m_plots.push_back( pair<string,string> ("nbjets", "# b jets") );
  clear();
}
//---------------------------------------------------------------------//
// Loop to make histograms
//---------------------------------------------------------------------//
void FakeClosurePlot::DataMCAnaPlots()
{
  // Here Loop over the channels and the signal regions that are set
  // via the run options.  Right now we loop over all and save all.
  string savedir(m_outputdir);
  for(uint ipr=0; ipr<m_PRs.size(); ++ipr){
    PlotRegion pr = m_PRs.at(ipr);
    string region = PRNames[pr];
    for(int ich=0; ich<Ch_N; ++ich){
      if(ich == Ch_all) continue;
      Chan ch        = (Chan) ich;
      string ch_name = chanNames[ch];
      for(uint ip=0; ip<m_plots.size(); ++ip){
        string var    = m_plots.at(ip).first;
        string xtitle = m_plots.at(ip).second;
        buildHists(m_hists, m_sys, var, xtitle, ch, pr);
        m_errs.push_back(buildErrors(m_hists.back(), m_sys));
        m_errs.push_back(buildRatioErrors(m_hists.back(), m_errs.back()));
        float xleg[] = {0.65,0.9};
        float yleg[] = {0.5,0.89};
        TLegend* leg = buildLegend(m_hists, m_errs[0], xleg, yleg);
        string save = savedir + region + "_" + ch_name + "_" + var + ".pdf";
        plotAll(m_hists, m_errs, save, leg, ich, false); //true);
        clear();
      }// end loop over plots
    }// end loop over channels
  }// end loop over regions
}

//---------------------------------------------------------------------//
// Set Histograms
//---------------------------------------------------------------------//
void FakeClosurePlot::buildHists(vector<TH1F*> &hists, vector<TH1F*> &sys, string var,
			       string xtitle, Chan ch, PlotRegion PR)
{
  if(m_dbg) cout << "FakeClosurePlot::buildHists from PR = " << PR << endl;
  string plot = PRNames[PR] + u() + chanNames[ch] + u() + var;
  if(m_dbg) cout << "Getting... " << plot <<endl;
  TH1F* SM = Get(m_files.at(0).file, plot + "_NOM", m_MCColor);
  SM->Reset();
  SM->SetName("SM");
  TH1F* sys_up = Get(m_files.at(0).file, plot + "_NOM", kBlack);
  sys_up->Reset();
  sys_up->SetName("sys_up");
  TH1F* sys_dn = Get(m_files.at(0).file, plot + "_NOM", kBlack);
  sys_dn->Reset();
  sys_dn->SetName("sys_dn");
  sys.push_back(sys_up);
  sys.push_back(sys_dn);
  // Loop over files and build histograms
  for(uint f=0; f<m_files.size(); ++f){
    File F = m_files.at(f);
    string end = !F.isfake ? "_NOM" : "_NONE";
    TH1F* h = Get(F.file, plot + end, F.color, xtitle.c_str(), "Entries");
    h->SetName(F.name.c_str());
    h->SetFillColor(F.color);
    hists.push_back(h);
    if(F.isfake)    addFakeSys(h, F.file, plot, sys); // Get Sys
    else if(F.ismc) addSysError(h, F.file, plot, sys);
    if( F.ismc ) SM->Add(h, 1.); // Add to SM
  }
  hists.push_back(SM);
}
//---------------------------------------------------------------------//
// Build legend
//---------------------------------------------------------------------//
TLegend* FakeClosurePlot::buildLegend(vector<TH1F*> hists, TGraphAsymmErrors* errs,
				    float* x, float* y)
{
  // All histograms that are in the list get put into the legend. Will
  // have to add the TGraphAsymmErrors later...
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
  return legend;

}
//---------------------------------------------------------------------//
// Build Ratio
//---------------------------------------------------------------------//
TH1F* FakeClosurePlot::buildRatio(TH1F* data, TH1F* SM)
{
  if(m_dbg) cout << "FakeClosurePlot::buildRatio" << endl;
  TH1F* ratio = RatioHist(data, SM, "Data/SM");
  return ratio;
}

//---------------------------------------------------------------------//
// Build the THStack from a set of histograms
//---------------------------------------------------------------------//
THStack* FakeClosurePlot::buildStack(vector<TH1F*> hists)
{
  if(m_dbg) cout << "FakeClosurePlot::buildStack" << endl;
  uint begin = hists.size() - 2; // Don't include total MC in the stack
  uint end   = 1;                // Don't include data in the stack
  THStack* stack = new THStack("stack","stack");
  for(uint i=begin; i>=end; i--) stack->Add( hists.at(i) );
  return stack;
}
//---------------------------------------------------------------------//
// Create and add error bars
//---------------------------------------------------------------------//
TGraphAsymmErrors* FakeClosurePlot::buildErrors(TH1F* summary, vector<TH1F*> sys)
{
  // Pass the summary histogram so we can get both the nominal values
  // and the statistical error from it. Also will add functions later
  // to incorporate the various sources of error.
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
TGraphAsymmErrors* FakeClosurePlot::buildRatioErrors(TH1F* nominal, TGraphAsymmErrors* tg_errs)
{
  // Pass the summary histogram so we can get both the nominal values
  // and the statistical error from it. Also will add functions later
  // to incorporate the various sources of error.
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
  }
  TGraphAsymmErrors* errs = new TGraphAsymmErrors(nbins, x, y, xerr, xerr, ym, yp);
  errs->SetMarkerSize(0);
  errs->SetFillStyle(3004);
  errs->SetFillColor(kGray+3);
  return errs;
}
//---------------------------------------------------------------------//
void FakeClosurePlot::addFakeSys(TH1F* nominal, TFile* file, string plot,
			       vector<TH1F*> &sys)
{
  // There are 8 total sys shifts that have the following names:
  // * {EL,MU}_RE_{UP,DOWN} = 4 shifts
  // * {EL,MU}_FR_{UP,DOWN} = 4 shifts
  // * Also statistical error

  vector<TH1F*> shifts; // Load Fake systematic shifts
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
void FakeClosurePlot::getFakeSys(TH1F* nominal, TFile* file, string plot,
			       float &sysup, float &sysdn)
{
  // There are 8 total sys shifts that have the following names:
  // * {EL,MU}_RE_{UP,DOWN} = 4 shifts
  // * {EL,MU}_FR_{UP,DOWN} = 4 shifts
  // * Also statistical error

  vector<TH1F*> shifts; // Load Fake systematic shifts
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
    float bc = nominal->GetBinContent(bin);
    for(uint s=0; s<shifts.size(); ++s){
      float shift = shifts.at(s)->GetBinContent(bin);
      if(shift > bc) sysup += pow(shift-bc,2);
      else           sysdn += pow(shift-bc,2);
    }// end loop over sys
  }// end loop over bins
  sysup = sqrt(sysup);
  sysdn = sqrt(sysdn);
  for(uint i=0; i<shifts.size(); ++i) shifts.at(i)->Delete();
  shifts.clear();
}
//---------------------------------------------------------------------//
void FakeClosurePlot::addSysError(TH1F* nominal, TFile* file, string plot,
				vector<TH1F*> &sys)
{
  // Will add the mc systematics later.  Right now I am only adding
  // the statistical error
  int nbins = nominal->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    float stat = pow(nominal->GetBinError(bin),2);
    float err_up = stat;
    float err_dn = stat;
    sys[0]->SetBinContent(bin, sys[0]->GetBinContent(bin) + err_up);
    sys[1]->SetBinContent(bin, sys[1]->GetBinContent(bin) + err_dn);
  }// end loop over bins
}
//---------------------------------------------------------------------//
// General tools
//---------------------------------------------------------------------//
void FakeClosurePlot::plotAll(vector<TH1F*> hists, vector<TGraphAsymmErrors*> errs,
			    string save, TLegend* leg, int ch, bool logy, bool logx)
{
  if(m_dbg) cout << "plotAll" << endl;
  TCanvas* c = makeCanvas("c");
  TPad* _pTop;
  TPad* _pBot;
  makePads(c,_pTop,_pBot);
  _pTop->SetBottomMargin(0.02);
  _pBot->SetTopMargin(0.015);
  if(logx){ _pTop->SetLogx(); _pBot->SetLogx(); }
  TH1F* data     = hists.at(0);   // Plotting three basic objects:
  THStack* stack = buildStack(hists);
  TH1F* SM       = hists.at(hists.size()-1);
  TH1F* ratio    = buildRatio(data, SM);
  TLine* nom = getLine(data, 1.00, kBlack, 1); // Create lines for ratio plot
  TLine* up50 = getLine(data, 1.5, kBlack, 2);
  TLine* dn50 = getLine(data, 0.5, kBlack, 2);
  c->cd();       // Set the attributes for the top Pad via the data histogram
  _pTop->Draw();
  _pTop->cd();
  stack->Draw("hist");
  stack->GetYaxis()->SetTitleSize(0.065);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.06);
  stack->GetXaxis()->SetLabelSize(0.0);
  stack->GetYaxis()->SetLabelFont(42);
  vector<float> minmax =  getMinMax(hists);
  if(!logy){ stack->SetMinimum(minmax[0]); stack->SetMaximum(minmax[1]); }
  else if(logy) { _pTop->SetLogy(); stack->SetMinimum(1e-4); stack->SetMaximum(minmax[1] * 100); }

  stack->Draw("hist"); // Draw the top canvas objects
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

  ratio->SetMaximum(2.0);  // Fix ratio plots to show +/- 100%
  ratio->SetMinimum(0.0);
  ratio->SetLabelFont(42,"X");
  ratio->SetLabelFont(42,"Y");
  ratio->GetYaxis()->CenterTitle();
  ratio->GetXaxis()->SetTitleOffset(1.0);
  ratio->GetXaxis()->SetTitleSize(0.12);

  c->cd();  // Draw the ratio plot on the bottom Will need to add the error object
  _pBot->Draw();
  _pBot->cd();
  ratio->Draw("ep");
  errs[1]->Draw("E2 same");
  nom->Draw("same");
  up50->Draw("same");
  dn50->Draw("same");
  ratio->Draw("ep same");
  c->SetFillColor(kWhite);
  c->Update();
  c->SaveAs( save.c_str() );
  delete c;
  delete lat;
}
//----------------------------------------------------//
float FakeClosurePlot::getNorm(TH1* h)
{
  return h->Integral(0,-1);
}
//-----------------------------------------------------//
float FakeClosurePlot::getStat(TH1* h, float low, float high)
{
  float be = 0;
  int min = h->FindBin(low);
  int max = h->FindBin(high);
  for(int bin=min; bin<=max; ++bin){ be += pow(h->GetBinError(bin),2); }
  return sqrt(be);
}
//-----------------------------------------------------//
void FakeClosurePlot::setMinMax(TH1* &h, float min, float max)
{
  h->SetMinimum(min);
  h->SetMaximum(max);
}
//-----------------------------------------------------//
TLine* FakeClosurePlot::getLine(TH1* h, float y, int color, int style)
{
  float x0     = h->GetBinCenter(1) - h->GetBinWidth(1)/2.;
  int finalBin = h->GetNbinsX();
  float x1     = h->GetBinCenter(finalBin) + h->GetBinWidth(finalBin)/2.;
  TLine* line  = makeLine(x0, x1, y, y, color, style);
  return line;
}
//-----------------------------------------------------//
float FakeClosurePlot::getMax(TH1F* h[], int n)
{
  float max = -999;
  for(int i=0; i<n; ++i){ if( max < h[i]->GetMaximum() ) max = h[i]->GetMaximum(); }
  return max;
}
//----------------------------------------------------------
FakeClosurePlot& FakeClosurePlot::setTag(const std::string &name)
{
  m_tag = name;
  return *this;
}
//----------------------------------------------------------
FakeClosurePlot& FakeClosurePlot::setInputDir(const std::string &dir)
{
  if(!dirExists(dir)) cout<<"Warning, invalid input dir '"<<dir<<"'"<<endl;
  m_inputdir = dir;
  return *this;
}
//----------------------------------------------------------
FakeClosurePlot& FakeClosurePlot::setOuputDir(const std::string &dir)
{
  if(dirExists(dir)) m_outputdir = dir;
  else {
    const bool dirWasCreated(mkdirIfNeeded(dir).size()>0);
    if(dirWasCreated) m_outputdir = dir;
    else              m_outputdir = "./";
  }
  return *this;
}
//----------------------------------------------------------
