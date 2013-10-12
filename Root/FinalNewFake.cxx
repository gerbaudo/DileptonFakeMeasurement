// This script will format the fake inputs to be used in
// the matrix method.

#include "SusyTest0/FinalNewFake.h"
#include "SusyTest0/utils.h"

#include<numeric> // accumulate

// same old ugly hack, to be removed ASAP. (DG 2013-09-17)
string getSrName(int sr) { return (SR_WHSS!=sr ? SRNames[sr] : CRNames[CR_SRWHSS]); }
string samplesHeader("\t\tttbar\t Wjet\t Zjet\t dib\t bbbar");
//----------------------------------------------------------
FinalNewFake::FinalNewFake(string outputname) :
  FakePlotting(RO_N),
  m_outfile(NULL)
{
}
//----------------------------------------------------------
void FinalNewFake::initIoFiles()
{
  string tag = "_"+m_tag;
  string inDir = m_inputdir+"/";
  initInputFile(m_mc,      inDir+"allBkg"     +tag+".root", "Combined MC", "mc",      kBlack,   25);
  initInputFile(m_ttbar,   inDir+"ttbar"      +tag+".root", "Ttbar",       "ttbar",   kBlue,    25);
  initInputFile(m_Wjet,    inDir+"wjets"      +tag+".root", "W+jet",       "wjet",    kMagenta, 25);
  initInputFile(m_Zjet,    inDir+"zjets"      +tag+".root", "Z+jet",       "Zjet",    kRed,     25);
  initInputFile(m_diboson, inDir+"diboson"    +tag+".root", "Diboson",     "diboson", kOrange,  22);
  initInputFile(m_bbbar,   inDir+"heavyflavor"+tag+".root", "b-bbar",      "bbbar",   kBlue,    23);
  m_files[FP_ttbar] = m_ttbar;
  m_files[FP_Wjet]  = m_Wjet;
  m_files[FP_Zjet]  = m_Zjet;
  m_files[FP_dib]   = m_diboson;
  m_files[FP_bbbar] = m_bbbar;
  m_outfile = new TFile(m_outputfname.c_str(), "recreate");
  cout<<"FinalNewFake: saving output to: '"<<m_outfile->GetName()<<"'"<<endl;
}
//----------------------------------------------------------
FinalNewFake::~FinalNewFake()
{
  m_mc.file->Close();
  m_ttbar.file->Close();
  m_Wjet.file->Close();
  m_Zjet.file->Close();
  m_diboson.file->Close();
}
//----------------------------------------------------------
void FinalNewFake::buildRates()
{
  buildMuonRateSR();
  buildElectronRateSR();
  buildSystematics();
  writeFileAndClose();
}
//----------------------------------------------------------
void clearTh1Vec(const vector<TH1*> &v)
{
  for(size_t i=0; i<v.size(); ++i) if(TH1 *h = v.at(i)) h->Delete();
}
//----------------------------------------------------------
string combinedHistoName(SignalRegion region, bool isEle, bool fake)
{
  string reg(SRNames[region]);
  string lep(isEle ? "el":"mu");
  string var(fake ? "fake_rate":"real_eff");
  return lep+"_"+var+"_"+reg;
}
//----------------------------------------------------------
string combinedHistoTitle(SignalRegion region, bool isEle, bool fake)
{
  string reg(SRProperNames[region]);
  string lep(isEle ? "Electron":"Muon");
  string var(fake ? "Fake Rate":"Real Eff");
  return lep+" "+var+": "+reg;
}
//----------------------------------------------------------
void printPercentages(const FinalNewFake::rlf_t config,
                      const vector<double> &v, const string &lab)
{
  cout<<(config.isFake?"Fake rate":"Real eff")
      <<": percentages for sr '"<<SRNames[config.reg]<<"'"
      <<", "<<(config.isEle?"electron":"muon")<<endl
      <<samplesHeader<<endl
      <<(config.isEle?"el":"mu")<<"_percent_"<<lab<<": "<<vdouble2str(v)<<endl
      <<endl;
}
//----------------------------------------------------------
void printPercentages(const FinalNewFake::rlf_t config,
                      const vector<double> &v1, const string &lab1,
                      const vector<double> &v2, const string &lab2)
{
  cout<<(config.isFake?"Fake rate":"Real eff")
      <<": percentages for sr '"<<SRNames[config.reg]<<"'"
      <<", "<<(config.isEle?"electron":"muon")<<endl
      <<samplesHeader<<endl
      <<(config.isEle?"el":"mu")<<"_percent_"<<lab1<<": "<<vdouble2str(v1)<<endl
      <<(config.isEle?"el":"mu")<<"_percent_"<<lab2<<": "<<vdouble2str(v2)<<endl
      <<endl;
}
//----------------------------------------------------------
bool FinalNewFake::combineWriteAndPlot(rlf_t cfg,
                                       cvecth1p_t &histos1, cvecth1p_t &histos2,
                                       cvecd_t &weights1, cvecd_t &weights2)
{
  if(TH1 *combinedHisto = getFinalRate(histos1, histos2, weights1, weights2)) {
    combinedHisto->SetName(combinedHistoName(cfg.reg, cfg.isEle, cfg.isFake).c_str());
    combinedHisto->SetTitle(combinedHistoTitle(cfg.reg, cfg.isEle, cfg.isFake).c_str());
    writeToOutputFile(combinedHisto);
    dumpPlot(combinedHisto);
    combinedHisto->Delete();
    return true;
  } else {
    cout<<"combineWriteAndPlot : invalid combinedHisto for "<<cfg.str()<<endl;
    return false;
  }
}
//----------------------------------------------------------
bool FinalNewFake::combineWriteAndPlot(rlf_t cfg, cvecth1p_t &histos, cvecd_t &weights)
{

  if(TH1 *combinedHisto = getFinalRate(histos, weights)) {
    combinedHisto->SetName(combinedHistoName(cfg.reg, cfg.isEle, cfg.isFake).c_str());
    combinedHisto->SetTitle(combinedHistoTitle(cfg.reg, cfg.isEle, cfg.isFake).c_str());
    writeToOutputFile(combinedHisto);
    dumpPlot(combinedHisto);
    combinedHisto->Delete();
    return true;
  } else {
    cout<<"combineWriteAndPlot : invalid combinedHisto for "<<cfg.str()<<endl;
    return false;
  }
}
//----------------------------------------------------------
void FinalNewFake::buildElectronRateSR()
{
  // Take into account the conv vs qcd in each signal region
  // so every signal region gets its own weighted fake rate
  const bool isEle(true);
  const string lepton("elec");
  vector<TH1*>   el_contrib_real, el_contrib_qcd, el_contrib_conv;
  for(int ic = 0; ic < FP_N; ++ic){ // Load the various rates from MC
    el_contrib_qcd .push_back(getFakeRate(lepton, (FakeProcess)ic, false));
    el_contrib_conv.push_back(getFakeRate(lepton, (FakeProcess)ic, true ));
    el_contrib_real.push_back(getRealEff (lepton, (FakeProcess)ic       ));
  }
  for(int sr = 0; sr<SR_N; ++sr){ // Fake Rate
    const bool isFake(true);
    string srname(getSrName(sr));
    vector<double> el_percent_qcd, el_percent_conv;
    rlf_t cfg(SignalRegion(sr), isEle, isFake);
    getFakePercentages(lepton, el_percent_qcd, el_percent_conv,  srname);
    printPercentages(cfg, el_percent_qcd, "qcd", el_percent_conv, "conv");
    if(combineWriteAndPlot(cfg, el_contrib_qcd, el_contrib_conv, el_percent_qcd, el_percent_conv)){
      if(m_dbg) cout<<"FinalNewFake::buildElectronRateSR : done fake for "<<cfg.str()<<endl;
    } else {
      cout<<"FinalNewFake::buildElectronRateSR : failed fake for "<<cfg.str()<<endl;
    }
  } // end for(sr)
  for(int sr = 0; sr<SR_N; ++sr){ // Real Efficiency
    const bool isFake(false);
    string srname(getSrName(sr));
    vector<double> el_percent_real;
    rlf_t cfg(SignalRegion(sr), isEle, isFake);
    getRealPercentages(lepton, el_percent_real,  srname);
    printPercentages(cfg, el_percent_real, "real");
    if(combineWriteAndPlot(cfg, el_contrib_real,el_percent_real)) {
      if(m_dbg) cout<<"FinalNewFake::buildElectronRateSR : done real for "<<cfg.str()<<endl;
    } else {
      cout<<"FinalNewFake::buildElectronRateSR : failed real for "<<cfg.str()<<endl;
    }
  } // end for(sr)
  clearTh1Vec(el_contrib_qcd);
  clearTh1Vec(el_contrib_conv);
  clearTh1Vec(el_contrib_real);
}
//----------------------------------------------------------
void FinalNewFake::buildMuonRateSR()
{
  // Take into account the conv vs qcd in each signal region
  // so every signal region gets its own weighted fake rate
  const bool isEle(false);
  string lepton = "muon";
  vector<TH1*>   mu_contrib_real, mu_contrib_qcd;
  for(int ic = 0; ic < FP_N; ++ic){ // Load the various rates from MC
    mu_contrib_qcd.push_back ( getFakeRate(lepton, (FakeProcess)ic, false) );
    mu_contrib_real.push_back( getRealEff (lepton, (FakeProcess)ic       ) );
  }
  for(int sr = 0; sr<SR_N; ++sr){ // Fake Rate
    const bool isFake(true);
    string srname(getSrName(sr));
    rlf_t cfg(SignalRegion(sr), isEle, isFake);
    vector<double> mu_percent_qcd, mu_percent_conv;
    getFakePercentages(lepton, mu_percent_qcd, mu_percent_conv,  srname);
    printPercentages(cfg, mu_percent_qcd, "qcd");
    if(combineWriteAndPlot(cfg, mu_contrib_qcd, mu_percent_qcd)) {
      if(m_dbg) cout<<"FinalNewFake::buildMuonRateSR : done fake for "<<cfg.str()<<endl;
    } else {
      cout<<"FinalNewFake::buildMuonRateSR : failed fake for "<<cfg.str()<<endl;
    }
  } // end for(sr)
  for(int sr = 0; sr<SR_N; ++sr){ // Real Efficiency
    const bool isFake(false);
    string srname(getSrName(sr));
    rlf_t cfg(SignalRegion(sr), isEle, isFake);
    vector<double> mu_percent_real;
    getRealPercentages(lepton, mu_percent_real,  srname);
    printPercentages(cfg, mu_percent_real, "real");
    if(combineWriteAndPlot(cfg, mu_contrib_real,mu_percent_real)) {
      if(m_dbg) cout<<"FinalNewFake::buildMuonRateSR : done real for "<<cfg.str()<<endl;
    } else {
      cout<<"FinalNewFake::buildMuonRateSR : failed real for "<<cfg.str()<<endl;
    }
  } // end for(sr)
  clearTh1Vec(mu_contrib_qcd);
  clearTh1Vec(mu_contrib_real);
}
//----------------------------------------------------------
void FinalNewFake::buildSystematics()
{
  // Loop over electron and muon systematics. In the end the Matrix method
  // code will combine them, but it is good in case we want to check
  // the effects of the individual systematics on the signal region
  vector<string> leptons;
  leptons.push_back("elec");
  leptons.push_back("muon");
  for(uint i=0; i<leptons.size(); ++i){
    string lep = leptons.at(i);
    pair< TParameter<double>, TParameter<double> > real = getRealSys(lep);
    //TH1*  metrel = getMetRelSys(lep);
    TH1* eta = getEtaSys(lep);
    TParameter<double> HFLFerr = getHFLFSys(lep);
    TParameter<double> datamc  = getDataMCSys(lep);
    TParameter<double> region  = getRegionSys(lep);
    //save(metrel);
    writeToOutputFile(eta);
    writeToOutputFile(real.first);
    writeToOutputFile(real.second);
    writeToOutputFile(HFLFerr);
    writeToOutputFile(datamc);
    writeToOutputFile(region);
    //metrel->Delete();
    eta->Delete();
  }
}
//----------------------------------------------------------
TH1* FinalNewFake::getRealEff(string lep, FakeProcess process)
{
  // Taken from MC in the data-driven Z control region.
  // Update: Just take data rate since we are scaling MC to data
  // and we are now wanting to apply Pt dependent sf to electrons.
  TH1 *rate = 0;
  string plot   = lep + "_realMC_all_l_pt_coarse";
  string ylabel = "Real Eff";
  File &f = m_files[process];
  if(f.file){
    rate    = buildRate(f.file, plot, f.sname, lep + " P_{T}", ylabel, kBlack);
    if( lep == "muon" ) scale(rate, mu_realSF);
    else                scale(rate, el_realSF);
  } else {
    cout<<"FinalNewFake::getRealEff : invalid file for '"<<process<<"'"<<endl;
  }
  return rate;
}
//----------------------------------------------------------
TH1* FinalNewFake::getFakeRate(string lep, FakeProcess process, bool isConv)
{
  // Now specify the file and the label
  string plot = lep + (isConv ? "_convMC_all_l_pt_coarse" : "_qcdMC_all_l_pt_coarse");
  //if( lep == "elec" && !isConv ) plot = lep + "_qcdMC_all_l_pt_eta";
  string ylabel = "Fake Rate";
  TH1* rate    = buildRate(m_files[process].file, plot,
                           m_files[process].sname,
                           lep + " P_{T}", ylabel, kBlack);
  if(rate){
    if(isConv) scale(rate, el_convSF);
    else if(lep == "elec") scale(rate, el_qcdSF);
    else if(lep == "muon") scale(rate, mu_qcdSF);
  } else {
    cout<<"FinalNewFake::getFakeRate: invalid rate histo "<<endl;
  }
  return rate;
}
//----------------------------------------------------------
void FinalNewFake::getFakePercentages(string lep,
                                      vector<double> &frac_qcd,
                                      vector<double> &frac_conv,
                                      string cr)
{
  // Loop over the regions and construct the percentages.
  // Currently we are taking everything for 40 < metrel < 100
  // from the OS Jet veto region, which is signal like. Percentages
  // are relative to the loose sample.
  // Currently only setup for electrons, but could be modified.
  frac_qcd.clear();
  frac_conv.clear();
  // Load the Loose Sample Numbers
  float total  = 0; // only one total for the final combination
  string plot = lep + "_" + cr + "_all_flavor_den";
  vector<double> temp_qcd;
  vector<double> temp_conv;
  for(int fp=0; fp<FP_N; ++fp){
    File &f = m_files[fp];
    if(TH1* h = Get(f.file, plot)) {
      float qcd(h->GetBinContent(LS_QCD+1)), conv(h->GetBinContent(LS_Conv+1));
      temp_qcd .push_back(qcd);
      temp_conv.push_back(conv);
      total += (qcd+ conv);
      h->Delete();
    } else {
      cout<<"FinalNewFake::getFakePercentages: cannot get '"<<plot<<"'"<<endl;
      temp_qcd .push_back(0.0);
      temp_conv.push_back(0.0);
    }
  }
  // -- DG you are here
  // -- float total(std::accumulate(temp_qcd.begin(), temp_qcd.end(), 0.0)
  // --             + std::accumulate(temp_conv.begin(), temp_conv.end(), 0.0));
  if(total == 0) cout << "In region: " << cr << " don't have any fakes!!! " << endl;
  // Set the percentage
  for(uint i=0; i<temp_qcd.size();  ++i) frac_qcd.push_back (total>0 ? temp_qcd .at(i)/total:0.);
  for(uint i=0; i<temp_conv.size(); ++i) frac_conv.push_back(total>0 ? temp_conv.at(i)/total:0.);
}
//----------------------------------------------------------
void FinalNewFake::getRealPercentages(string lep,
                                      vector<double> &frac,
                                      string cr)
{
  // Loop over the regions and construct the percentages.
  // Currently we are taking everything for 40 < metrel < 100
  // from the OS Jet veto region, which is signal like. Percentages
  // are relative to the loose sample.
  // Currently only setup for electrons, but could be modified.
  frac.clear();
  float total  = 0; // only one total for the final combination
  vector<double> temp;
  string plot  = lep + "_" + cr + "_all_flavor_den";
  for(int fp=0; fp<FP_N; ++fp){
    if(TH1* h = Get(m_files[fp].file, plot)) {
      float contrib  = h->GetBinContent(LS_Real+1);
      temp.push_back( contrib );
      total += contrib;
      h->Delete();
    } else {
      cout<<"FinalNewFake::getRealPercentages: cannot get '"<<plot<<"'"<<endl;
      temp.push_back(0.0);
    }
  } // end for(fp)
  // -- DG you are here
  // -- float total(std::accumulate(temp.begin(), temp.end(), 0.0));
  if(total == 0) cout << "In region: " << cr << " don't have any reals!!! " << endl;
  for(uint i = 0; i<temp.size(); ++i) frac.push_back( total > 0 ? temp.at(i) / total : 0. );
}
//----------------------------------------------------------
TH1* cloneLast(const vector<TH1*> hs, const string &name)
{
  TH1 *h = NULL;
  if(hs.back()) h = static_cast<TH1*>(hs.back()->Clone(name.c_str()));
  else cout<<"cloneLast: invalid last histo"<<endl;
  return h;
}
//----------------------------------------------------------
struct SumWithErr2 {
  float sum, err2;
  SumWithErr2(float s, float e2) : sum(s), err2(e2) {}
};
//----------------------------------------------------------
SumWithErr2 binWeightedSum(const vector<TH1*> &hs, const vector<double> &ws, int bx)
{
  SumWithErr2 se2(0.0, 0.0);
  if(hs.size()!=ws.size()) cout<<"binWeightedSum invalid inputs"<<hs.size()<<","<<ws.size()<<endl;
  else {
    for(uint i=0; i<hs.size(); ++i){
      if(TH1 *h = hs.at(i)){
        float bc(h->GetBinContent(bx)), be(h->GetBinError(bx)), w(ws.at(i));
        se2.sum  += bc * w;
        se2.err2 += be * be * w * w;
      } else { cout<<"binWeightedSum invalid h["<<i<<"]"<<endl; continue; }
    }
  }
  return se2;
}
//----------------------------------------------------------
// DG : this was probably used for some 2d-parametrization test; currently not needed
SumWithErr2 bin2dWeightedSum(const vector<TH1*> &hs, const vector<double> &ws, int bx, int by)
{
  SumWithErr2 se2(0.0, 0.0);
  if(hs.size()!=ws.size()) cout<<"bin2dWeightedSum invalid inputs"<<hs.size()<<","<<ws.size()<<endl;
  else {
    for(uint i=0; i<hs.size(); ++i){
      if(TH1 *h = hs.at(i)){
        float bc(h->GetBinContent(bx,by)), be(h->GetBinError(bx, by)), w(ws.at(i));
        se2.sum  += bc * w;
        se2.err2 += be * be * w * w;
      } else { cout<<"bin2dWeightedSum invalid h["<<i<<"]"<<endl; continue; }
    }
  }
  return se2;
}
//----------------------------------------------------------
TH1* FinalNewFake::getFinalRate(vector<TH1*> rates_qcd,
                                vector<TH1*> rates_conv,
                                vector<double> percent_qcd,
                                vector<double> percent_conv)
{
  // Here combine the rates based on their relative percentages.
  // It is assumed that the ith component of rates corresponds
  // to the ith component of percentages.
  TH1* final = NULL;
  if((final = cloneLast(rates_qcd, "final_rate"))) {
    final->Reset();
    int nbinsx(final->GetNbinsX());
    for(int bx=1; bx<=nbinsx; ++bx){
      SumWithErr2 qcd (binWeightedSum(rates_qcd,  percent_qcd,  bx));
      SumWithErr2 conv(binWeightedSum(rates_conv, percent_conv, bx));
      final->SetBinContent(bx, qcd.sum+conv.sum );
      final->SetBinError(bx, TMath::Sqrt(qcd.err2+conv.err2));
    }// end for(bx)
  } else {
    cout<<"FinalNewFake::getFinalRate: cannot get template histo"<<endl;
  }
  return final;
}
//----------------------------------------------------------
TH1* FinalNewFake::getFinalRate(vector<TH1*> rates,
                                vector<double> percent)
{
  // Here combine the rates based on their relative percentages.
  // It is assumed that the ith component of rates corresponds
  // to the ith component of percentages.
  TH1* final = NULL;
  if((final = cloneLast(rates, "final_rate"))) {
    final->Reset();
    int nbins = final->GetNbinsX();
    for(int bx=1; bx<=nbins; ++bx){ // Loop over the bins and weight
      SumWithErr2 tot(binWeightedSum(rates, percent, bx));
      final->SetBinContent(bx, tot.sum );
      final->SetBinError(bx, TMath::Sqrt(tot.err2));
    } // end for(bx)
  } else {
    cout<<"FinalNewFake::getFinalRate: cannot clone '"<<rates.back()->GetName()<<"'"<<endl;
  }
  return final;
}
//----------------------------------------------------------
void FinalNewFake::scale(TH1* &h, float sf)
{
  // Scale each bin by the scale factor
  for(int bin=1; bin<=h->GetNbinsX(); ++bin){
    float bc = h->GetBinContent(bin);
    h->SetBinContent(bin, bc * sf );
  }
}
//----------------------------------------------------------
pair< TParameter<double>,TParameter<double> > FinalNewFake::getRealSys(string lep)
{
  // This is the combined real systeamtic.  Since it is small, not going
  // to divide it up into separate factors
  bool isElec(lep=="elec");
  string name(isElec ? "el" : "mu");
  TParameter<double> shiftup  ((name + "_real_up"  ).c_str(), isElec ? el_real_up : mu_real_up);
  TParameter<double> shiftdown((name + "_real_down").c_str(), isElec ? el_real_dn : mu_real_dn);
  return pair< TParameter<double>,TParameter<double> > (shiftup, shiftdown);
}
//----------------------------------------------------------
TParameter<double> FinalNewFake::getHFLFSys(string lep)
{
  bool isElec(lep == "elec");
  string name(isElec ? "el" : "mu");
  TParameter<double> shift((name + "_HFLFerr").c_str(), isElec ? el_HFLFerr : mu_HFLFerr);
  return shift;
}
//----------------------------------------------------------
TParameter<double> FinalNewFake::getDataMCSys(string lep)
{
  bool isElec(lep == "elec");
  string name(isElec ? "el" : "mu");
  TParameter<double> shift((name + "_datamc").c_str(), isElec ? el_datamc : mu_datamc);
  return shift;
}
//----------------------------------------------------------
TParameter<double> FinalNewFake::getRegionSys(string lep)
{
  bool isElec(lep == "elec");
  string name(isElec ? "el" : "mu");
  TParameter<double> shift((name+"_region").c_str(), isElec ? el_region : mu_region);
  return shift;
}
//----------------------------------------------------------
TH1* FinalNewFake::getMetRelSys(string lep)
{
  // Here we take the MetRel distribution (no met cut) and
  // normalize this to the average fake rate
  // (taken from one bin rate) and then take the differences from
  // 1 as the percentage.
  bool isElec(lep == "elec");
  string l(isElec ? "el" : "mu");
  string sn(m_mc.sname);
  TFile *f = m_mc.file;
  TH1* rate = buildRate(f, lep+"_qcdMC_all_metrel", sn, "#slash{E}^{rel}_{T}", "", kBlack);
  TH1* norm = buildRate(f, lep+"_qcdMC_all_onebin", sn, "onebin",              "", kBlack);
  rate->Scale(1/norm->GetBinContent(1));
  norm->Delete();
  rate->SetName((l+"_metrel_sys").c_str());
  rate->SetTitle("#slash{E}^{rel}_{T} systematic");
  for(int bin = 1; bin<= rate->GetNbinsX(); ++bin) rate->AddBinContent(bin, -1.0);
  return rate;
}
//----------------------------------------------------------
TH1* FinalNewFake::getEtaSys(string lep)
{
  // Here we take the eta distribution and
  // normalize this to the average fake rate
  // (taken from one bin rate) and then take the differences from
  // 1 as the percentage.
  bool isElec(lep == "elec");
  string name(isElec ? "el" : "mu");
  string outname(name+"_eta_sys");
  string sn(m_mc.sname);
  TFile *f = m_mc.file;
  TH1* rate = buildRate(f, lep+"_qcdMC_all_l_eta_coarse",  sn, "|#eta|", "", kBlack);
  TH1* norm = buildRate(f, lep+"_qcdMC_all_onebin",        sn, "onebin", "", kBlack);
  rate->Scale(1/norm->GetBinContent(1));
  rate->SetName(outname.c_str());
  norm->Delete();
  for(int bin=1; bin<=rate->GetNbinsX(); ++bin)
    if(isElec) rate->AddBinContent(bin, -1.0);
    else       rate->SetBinContent(bin,  0.0); // mu consistent with 0
  return rate;
}
//----------------------------------------------------------
FinalNewFake& FinalNewFake::setTag(const std::string &name)
{
  m_tag = name;
  return *this;
}
//----------------------------------------------------------
FinalNewFake& FinalNewFake::setInputDir(const std::string &dir)
{
  if(!dirExists(dir)) cout<<"Warning, invalid input dir '"<<dir<<"'"<<endl;
  m_inputdir = dir;
  return *this;
}
//----------------------------------------------------------
FinalNewFake& FinalNewFake::setOuputFilename(const std::string &name)
{
  const string dir(basedir(name));
  if(dirExists(dir)) m_outputfname = name;
  else {
    const bool dirWasCreated(mkdirIfNeeded(dir).size()>0);
    if(dirWasCreated) m_outputfname = name;
    else {
      const string fallbackDir("./");
      m_outputfname = name;
      replace(m_outputfname, dir, fallbackDir);
      cout<<"output path '"<<name<<"', output stored to '"<<m_outputfname<<"'"<<endl;
    }
  }
  return *this;
}
//----------------------------------------------------------
FinalNewFake& FinalNewFake::setOuputPlotdir(const std::string &name)
{
  m_outputplotdir = mkdirIfNeeded(name);
  bool invalidDir(m_outputplotdir.size()==0);
  if(invalidDir) {
    cout<<"invalid outputplotdir '"<<name<<"' trying to guess something reasonable..."<<endl;
    m_outputplotdir = basedir(m_outputfname);
    invalidDir = m_outputplotdir.size()==0;
    if(invalidDir) m_outputplotdir = "./";
  }
  if(!endswith(m_outputplotdir,"/")) m_outputplotdir.append("/");
  return *this;
}
//----------------------------------------------------------
void FinalNewFake::writeToOutputFile(TH1* h)
{
  m_outfile->cd();
  h->Write();
}
//----------------------------------------------------------
void FinalNewFake::writeToOutputFile(TParameter<double> param)
{
  m_outfile->cd();
  param.Write();
}

//----------------------------------------------------------
void FinalNewFake::writeFileAndClose()
{
  m_outfile->cd();
  m_outfile->Write();
  m_outfile->Close();
}
//----------------------------------------------------------
//------------------------------------------------------------//
// Simple plotting method
//------------------------------------------------------------//
void FinalNewFake::dumpPlot(TH1* hist)
{
  TCanvas* c = makeCanvas("c");
  float xleg[] = {0.15, 0.3}, yleg[] = {0.85, 0.9};
  TLegend* leg = makeLegend(xleg,yleg);
  leg->SetTextSize(0.04);
  leg->AddEntry(hist);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleOffset(0.85);
  hist->SetLabelSize(0.04, "Y");
  hist->SetLabelSize(0.04, "X");
  int bin(hist->GetMaximumBin());
  float max(hist->GetBinContent(bin) + hist->GetBinError(bin));
  hist->SetMaximum( max * 1.1 );
  hist->Draw("p");
  leg->Draw("same");
  string plotname(m_outputplotdir+hist->GetName());
  c->SaveAs((plotname+".eps").c_str());
  c->SaveAs((plotname+".png").c_str());
  delete c;
}
//----------------------------------------------------------
void FinalNewFake::dumpStat(TH1* hist)
{
  cout<<"Name: "<<hist->GetName();
  float stat = 0;
  for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
    float dev = (hist->GetBinError(bin)/hist->GetBinContent(bin)) * 100;
    if(stat < dev) stat = dev;
  }
  cout<<" Maximum deviation: "<<stat<<endl;
}
//----------------------------------------------------------
