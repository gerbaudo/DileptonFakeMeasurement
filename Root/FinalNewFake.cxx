// This script will format the fake inputs to be used in 
// the matrix method.

#include "SusyTest0/FinalNewFake.h"
#include "SusyTest0/utils.h"

// same old ugly hack, to be removed ASAP. (DG 2013-09-17)
string getSrName(int sr) { return (SR_WHSS!=sr ? SRNames[sr] : CRNames[CR_SRWHSS]); }
string samplesHeader("\t\tttbar\t Wjet\t Zjet\t dib\t bbbar");
//------------------------------------------------------------//
// Constructor
//------------------------------------------------------------//
FinalNewFake::FinalNewFake(string outputname) :
  FakePlotting(RO_N),
  m_outfile(NULL)
{
}
//----------------------------------------------------------
void FinalNewFake::initIoFiles()
{
  string tag = m_tag;
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

//------------------------------------------------------------//
// Destructor
//------------------------------------------------------------//
FinalNewFake::~FinalNewFake()
{
  m_mc.file->Close();
  m_ttbar.file->Close();
  m_Wjet.file->Close();
  m_Zjet.file->Close();
  m_diboson.file->Close();
}

//------------------------------------------------------------//
// Build Rates
//------------------------------------------------------------//
void FinalNewFake::buildRates()
{

  buildMuonRateSR();
  buildElectronRateSR();
  buildSystematics();
  write();

}

//------------------------------------------------------------//
// Build Electron Rate using signal region composition
//------------------------------------------------------------//
void FinalNewFake::buildElectronRateSR()
{
  // Take into account the conv vs qcd in each signal region
  // so every signal region gets its own weighted fake rate
  string lepton = "elec";
  vector<TH1*> el_contrib_real;
  vector<double> el_percent_real;
  vector<TH1*> el_contrib_qcd;
  vector<double> el_percent_qcd;  
  vector<TH1*> el_contrib_conv;
  vector<double> el_percent_conv;  
  TH1* el_fr = NULL;
  TH1* el_re = NULL;  
  // Load the various rates from MC
  for(int ic = 0; ic < FP_N; ++ic){
    el_contrib_qcd.push_back( getFakeRate(lepton, (FakeProcess)ic, false) );
    el_contrib_conv.push_back( getFakeRate(lepton, (FakeProcess)ic, true) );
    el_contrib_real.push_back( getRealEff(lepton, (FakeProcess)ic) );
  }
  // Create and Save Fake Rate
  for(int sr = 0; sr<SR_N; ++sr){
    string srname = getSrName(sr);
    if(m_dbg) cout<<"Getting for sr: "<<srname<<endl;
    el_percent_qcd.clear();
    el_percent_conv.clear();    
    getPercentages(lepton, el_percent_qcd, el_percent_conv,  srname);
    cout<<"Fake rate: percentages for sr '"<<srname<<"', electron"<<endl
        <<samplesHeader<<endl
        <<"el_percent_qcd: "<<vdouble2str(el_percent_qcd)<<endl
        <<"el_percent_conv: "<<vdouble2str(el_percent_conv)<<endl
        <<endl;
    if((el_fr = getFinalRate(el_contrib_qcd, el_contrib_conv, el_percent_qcd, el_percent_conv))) {
      el_fr->SetName(("el_fake_rate_"+srname).c_str());
      el_fr->SetTitle(("Electron Fake Rate: " + SRProperNames[sr]).c_str());
      save(el_fr);
      dumpPlot(el_fr);
      el_fr->Delete();
      el_fr = NULL;
    } else {
     cout<<"FinalNewFake::buildElectronRateSR : invalid el_fr for '"<<srname<<"'"<<endl;
    }
  } // end for(sr)
  // Clean up the fake rates
  for(uint i=0; i<el_contrib_qcd.size(); ++i)  if(TH1 *h = el_contrib_qcd.at(i))  h->Delete();
  for(uint i=0; i<el_contrib_conv.size(); ++i) if(TH1 *h = el_contrib_conv.at(i)) h->Delete();
  // Real Efficiency
  for(int sr = 0; sr<SR_N; ++sr){
    string srname = getSrName(sr);
    el_percent_real.clear();
    getPercentages(lepton, el_percent_real,  srname);
    cout<<"Real eff: percentages for sr '"<<srname<<"', electron"<<endl
        <<samplesHeader<<endl
        <<"el_percent_real: "<<vdouble2str(el_percent_real)<<endl
        <<endl;
    if(m_dbg) cout<<"Getting for sr: "<<srname<<endl;
    if((el_re = getFinalRate(el_contrib_real,el_percent_real))) {
      el_re->SetName(("el_real_eff_"+srname).c_str());
      el_re->SetTitle(("Electron Real Eff: " + SRProperNames[sr]).c_str());
      save(el_re);
      dumpPlot(el_re);
      el_re->Delete();
      el_re = NULL;
    } else {
      cout<<"FinalNewFake::buildElectronRateSR : invalid el_re for '"<<srname<<"'"<<endl;
    }
  } // end for(sr)
  // Clean up the fake rates
  for(uint i=0; i<el_contrib_real.size(); ++i) if(TH1 *h = el_contrib_real.at(i)) h->Delete();

}
//------------------------------------------------------------//
// Build Muon SR Rate
//------------------------------------------------------------//
void FinalNewFake::buildMuonRateSR()
{

  // Take into account the conv vs qcd in each signal region
  // so every signal region gets its own weighted fake rate
  
  string lepton = "muon";
  vector<TH1*> mu_contrib_real;
  vector<double> mu_percent_real;
  vector<TH1*> mu_contrib_qcd;
  vector<double> mu_percent_qcd;
  vector<TH1*> mu_contrib_conv;  // Dummy holder 
  vector<double> mu_percent_conv; // Dummy holder
  TH1* mu_fr = NULL;
  TH1* mu_re = NULL; 
  // Load the various rates from MC
  for(int ic = 0; ic < FP_N; ++ic){
    mu_contrib_qcd.push_back( getFakeRate(lepton, (FakeProcess)ic, false) );
    mu_contrib_real.push_back( getRealEff(lepton, (FakeProcess)ic) );
  }
  // Create and Save Fake Rate
  for(int sr = 0; sr<SR_N; ++sr){
    string srname = getSrName(sr);
    if(m_dbg) cout<<"Getting for sr: "<<srname<<endl;
    mu_percent_qcd.clear();
    mu_percent_conv.clear();    
    getPercentages(lepton, mu_percent_qcd, mu_percent_conv,  srname);
    cout<<"Fake rate: percentages for sr '"<<srname<<"', muon"<<endl
        <<samplesHeader<<endl
        <<"mu_percent_qcd: "<<vdouble2str(mu_percent_qcd)<<endl
        <<"mu_percent_conv: "<<vdouble2str(mu_percent_conv)<<endl
        <<endl;
    mu_fr = getFinalRate(mu_contrib_qcd, mu_contrib_conv, mu_percent_qcd, mu_percent_conv);
    if(mu_fr) {
      mu_fr->SetName(("mu_fake_rate_"+srname).c_str());
      mu_fr->SetTitle(("Muon Fake Rate: " + SRProperNames[sr]).c_str());
      save(mu_fr);
      dumpPlot(mu_fr);
      mu_fr->Delete();
      mu_fr = NULL;
    } else {
      cout<<"FinalNewFake::buildMuonRateSR : invalid mu_fr for '"<<srname<<"'"<<endl;
    }
  } // end for(sr)
  // Clean up the fake rates
  for(uint i=0; i<mu_contrib_qcd.size(); ++i) if(TH1* h = mu_contrib_qcd.at(i)) h->Delete();
  mu_contrib_conv.clear();

  // Real Efficiency
  for(int sr = 0; sr<SR_N; ++sr){
    string srname = getSrName(sr);
    mu_percent_real.clear();
    getPercentages(lepton, mu_percent_real,  srname);
    cout<<"Real eff: percentages for sr '"<<srname<<"', muon"<<endl
        <<samplesHeader<<endl
        <<"mu_percent_real: "<<vdouble2str(mu_percent_real)<<endl
        <<endl;
    if(m_dbg) cout<<"Getting for sr: "<<srname<<endl;
    mu_re = getFinalRate(mu_contrib_real,mu_percent_real);
    if(mu_re) {
      mu_re->SetName(("mu_real_eff_"+srname).c_str());
      mu_re->SetTitle(("Muon Real Eff: " + SRProperNames[sr]).c_str());
      dumpPlot(mu_re);
      save(mu_re);
      mu_re->Delete();
      mu_re = NULL;
    } else {
      cout<<"FinalNewFake::buildMuonRateSR : invalid mu_re for '"<<srname<<"'"<<endl;      
    }
  } // end for(sr)
  // Clean up the fake rates
  for(uint i=0; i<mu_contrib_real.size(); ++i) if(TH1* h = mu_contrib_real.at(i)) h->Delete();
}

//------------------------------------------------------------//
// Build Systematics
//------------------------------------------------------------//
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
    save(eta);
    save(real.first);
    save(real.second);
    save(HFLFerr);
    save(datamc);
    save(region);

    //metrel->Delete();
    eta->Delete();
  }

}

//------------------------------------------------------------//
// Get Real Efficiency
//------------------------------------------------------------//
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

//------------------------------------------------------------//
// Get Fake Rate
//------------------------------------------------------------//
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

//------------------------------------------------------------//
// Get the relative percentages and store in a histo
//------------------------------------------------------------//
void FinalNewFake::getPercentages(string lep, 
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
  // TODO: Revisit these regions to loosen cuts so these
  // special cases do not have to be used
  //if(cr == "SRmT2b"){cout<<"Reseting "<<cr<<endl; cr = "SRmT2a";}
  //if(cr == "SRmT2c"){cout<<"Reseting "<<cr<<endl; cr = "SRmT2a";}

  // Load the Loose Sample Numbers
  float total  = 0; // only one total for the final combination
  string plot = lep + "_" + cr + "_all_flavor_den"; 
  vector<double> temp_qcd;  
  vector<double> temp_conv;
  for(int fp=0; fp<FP_N; ++fp){
    File &f = m_files[fp];
    if(TH1* h = Get(f.file, plot)) {
      float qcd  = h->GetBinContent(LS_QCD+1);
      float conv = h->GetBinContent(LS_Conv+1);
      temp_qcd.push_back( qcd );
      temp_conv.push_back( conv );
      total += qcd;
      total += conv;
      h->Delete();
    } else {
      cout<<"FinalNewFake::getPercentages: cannot get '"<<plot<<"'"<<endl;
      temp_qcd.push_back(0.0);
      temp_conv.push_back(0.0);      
    }
  }
  if(total == 0) cout << "In region: " << cr << " don't have any fakes!!! " << endl;

  // Set the Percentage histograms
  for(uint i = 0; i<temp_qcd.size(); ++i) frac_qcd.push_back(total>0 ? temp_qcd.at(i)/total:0.);
  for(uint i = 0; i<temp_conv.size(); ++i)frac_conv.push_back(total>0 ? temp_conv.at(i)/total:0.);
}

//------------------------------------------------------------//
void FinalNewFake::getPercentages(string lep, 
                                  vector<double> &frac,
                                  string cr)
{

  // Loop over the regions and construct the percentages. 
  // Currently we are taking everything for 40 < metrel < 100
  // from the OS Jet veto region, which is signal like. Percentages
  // are relative to the loose sample.
  // Currently only setup for electrons, but could be modified.

  frac.clear();

  // Load the Loose Sample Numbers
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
      cout<<"FinalNewFake::getPercentages(: cannot get '"<<plot<<"'"<<endl;
      temp.push_back(0.0);
    }
  } // end for(fp)
  if(total == 0) cout << "In region: " << cr << " don't have any reals!!! " << endl;
  for(uint i = 0; i<temp.size(); ++i) frac.push_back( total > 0 ? temp.at(i) / total : 0. );
}

//------------------------------------------------------------//
// Get Final weighted histogram
//------------------------------------------------------------//
TH1* FinalNewFake::getFinalRate(vector<TH1*> rates, vector<TH1*> percentages)
{
  // Here combine the rates based on their relative percentages.
  // It is assumed that the ith component of rates corresponds
  // to the ith component of percentages.
  TH1* final = 0;
  if(!rates.back()) { cout<<"FinalNewFake::getFinalRate: invalid last histo"<<endl; return final; }
  if((final = static_cast<TH1*>(rates.back()->Clone("final_rate")))) {
    final->Reset();
    int nbins = final->GetNbinsX();
    for(int bin=1; bin<=nbins; ++bin){ // Loop over the bins and weight
      int ntypes = 0;
      float rate = 0;
      float err  = 0;
      if(rates.size()!=percentages.size()) {
        cout<<"getFinalRate:"
            <<" rates["<<rates.size()<<"]!=percentages["<<percentages.size()<<"]...skip"<<endl;
        continue;
      }
      for(uint i=0; i<rates.size(); ++i){
        float bc  = rates.at(i)->GetBinContent(bin);
        float be  = rates.at(i)->GetBinError(bin); 
        float per = percentages.at(i)->GetBinContent(bin);
        rate += bc * per;
        err  += be * per * be * per;
        ntypes++;
      }
      final->SetBinContent(bin,  ntypes == 0 ? 0.0 : rate / ntypes );
      final->SetBinError(bin,  ntypes == 0 ? 0.0 : TMath::Sqrt(err) / ntypes );
    }// end for(bin)
  } else {
    cout<<"FinalNewFake::getFinalRate: cannot clone '"<<rates.back()->GetName()<<"'"<<endl;
  }
  return final;

}
//------------------------------------------------------------//
TH1* FinalNewFake::getFinalRate(vector<TH1*> rates_qcd,
                                vector<TH1*> rates_conv,
                                vector<double> percent_qcd,
                                vector<double> percent_conv)
{
  // Here combine the rates based on their relative percentages.
  // It is assumed that the ith component of rates corresponds
  // to the ith component of percentages.
  TH1* final = 0;
  if(!rates_qcd.back()) { cout<<"FinalNewFake::getFinalRate(: invalid last histo"<<endl; return final; }
  if((final = static_cast<TH1*>(rates_qcd.back()->Clone("final_rate")))) {
    final->Reset();
    int nbinsx(final->GetNbinsX()), nbinsy(final->GetYaxis()->GetNbins());
    for(int binx=1; binx<=nbinsx; ++binx){ // Loop over the bins and weight
      for(int biny=1; biny<=nbinsy; ++biny){
        float rate = 0;
        float err  = 0;
        for(uint i=0; i<rates_qcd.size(); ++i){ // QCD Rate
          TH1 *hqcd = rates_qcd.at(i);
          bool sameSize(rates_qcd.size()==percent_qcd.size());
          if(hqcd && sameSize){
            float bc(hqcd->GetBinContent(binx,biny)), be(hqcd->GetBinError(binx, biny)); 
            float per = percent_qcd.at(i);
            rate += bc * per;
            err  += be * per * be * per;
          } else {
            cout<<"invalid inputs "<<hqcd<<", "<<rates_qcd.size()<<"!="<<percent_qcd.size()<<endl;
            continue;
          }
        } // end for(i)
        for(uint i=0; i<rates_conv.size(); ++i){ // Conv Rate
          TH1 *hcon = rates_conv.at(i);
          bool sameSize(rates_conv.size()==percent_conv.size());
          if(hcon && sameSize){
            float bc(hcon->GetBinContent(binx, biny)), be(hcon->GetBinError(binx, biny));
            float per = percent_conv.at(i);
            rate += bc * per;
            err  += be * per * be * per;
          } else {
            cout<<"invalid inputs "<<hcon<<", "<<rates_conv.size()<<"!="<<percent_conv.size()<<endl;
            continue;
          }
        }
        final->SetBinContent(binx, biny, rate );
        final->SetBinError(binx, biny, TMath::Sqrt(err));
      }// end for(biny)
    }// end for(binx)
  } else {
    cout<<"FinalNewFake::getFinalRate: cannot clone '"<<rates_conv.back()->GetName()<<"'"<<endl;    
  }
  return final;
}
//------------------------------------------------------------//
TH1* FinalNewFake::getFinalRate(vector<TH1*> rates,
                                vector<double> percent)
{
  // Here combine the rates based on their relative percentages.
  // It is assumed that the ith component of rates corresponds
  // to the ith component of percentages.
  TH1* final = 0;
  if(!rates.back()) { cout<<"FinalNewFake::getFinalRate: invalid last histo"<<endl; return final; }
  if((final = static_cast<TH1*>(rates.back()->Clone("final_rate")))) {
    final->Reset();
    int nbins = final->GetNbinsX();
    for(int bin=1; bin<=nbins; ++bin){ // Loop over the bins and weight
      float rate = 0;
      float err  = 0;
      if(rates.size()!=percent.size()) {
        cout<<"getFinalRate:"
            <<" rates["<<rates.size()<<"]!=percent["<<percent.size()<<"]...skip"<<endl;
        continue;
      }
      for(uint i=0; i<rates.size(); ++i){ // Rate
        if(TH1* hr = rates.at(i)){
          float bc(hr->GetBinContent(bin)), be(hr->GetBinError(bin));
          float per = percent.at(i);
          rate += bc * per;
          err  += be * per * be * per;
        } else { cout<<"getFinalRate invalid TH1* rates["<<i<<"]"<<endl; continue; }
      } // end for(i)
      final->SetBinContent(bin, rate );
      final->SetBinError(bin, TMath::Sqrt(err));
    } // end for(bin)
  } else {
    cout<<"FinalNewFake::getFinalRate: cannot clone '"<<rates.back()->GetName()<<"'"<<endl;
  }
  return final;
}


//------------------------------------------------------------//
// Scale histogram by constant factor
//------------------------------------------------------------//
void FinalNewFake::scale(TH1* &h, float sf)
{

  // Scale each bin by the scale factor
  for(int bin=1; bin<=h->GetNbinsX(); ++bin){
    float bc = h->GetBinContent(bin);
    h->SetBinContent(bin, bc * sf );
  }

}

//------------------------------------------------------------//
// Get Systematic variations
//------------------------------------------------------------//
//------------------------------------------------------------//
pair< TParameter<double>,TParameter<double> > FinalNewFake::getRealSys(string lep)
{

  // This is the combined real systeamtic.  Since it is small, not going
  // to divide it up into separate factors

  bool isElec = lep == "elec";
  string name = isElec ? "el" : "mu";

  TParameter<double> shiftup = isElec ? 
    TParameter<double> ((name + "_real_up").c_str(),el_real_up) : 
    TParameter<double> ((name + "_real_up").c_str(),mu_real_up);

  TParameter<double> shiftdown = isElec ? 
    TParameter<double> ((name + "_real_down").c_str(),el_real_dn) : 
    TParameter<double> ((name + "_real_down").c_str(),mu_real_dn);

  return pair< TParameter<double>,TParameter<double> > (shiftup, shiftdown);

}

//------------------------------------------------------------//
TParameter<double> FinalNewFake::getHFLFSys(string lep)
{

  bool isElec = lep == "elec";
  string name = isElec ? "el" : "mu";
  
  TParameter<double> shift = isElec ?
    TParameter<double> ((name + "_HFLFerr").c_str(), el_HFLFerr) :
    TParameter<double> ((name + "_HFLFerr").c_str(), mu_HFLFerr);

  return shift;
}

//------------------------------------------------------------//
TParameter<double> FinalNewFake::getDataMCSys(string lep)
{

  bool isElec = lep == "elec";
  string name = isElec ? "el" : "mu";
  
  TParameter<double> shift = isElec ?
    TParameter<double> ((name + "_datamc").c_str(), el_datamc) :
    TParameter<double> ((name + "_datamc").c_str(), mu_datamc);

  return shift;
}
//------------------------------------------------------------//
TParameter<double> FinalNewFake::getRegionSys(string lep)
{

  bool isElec = lep == "elec";
  string name = isElec ? "el" : "mu";
  
  TParameter<double> shift = isElec ?
    TParameter<double> ((name + "_region").c_str(), el_region) :
    TParameter<double> ((name + "_region").c_str(), mu_region);

  return shift;
}
//------------------------------------------------------------//
TH1* FinalNewFake::getMetRelSys(string lep)
{
  // Here we take the MetRel distribution (no met cut) and 
  // normalize this to the average fake rate 
  // (taken from one bin rate) and then take the differences from
  // 1 as the percentage.

  bool isElec = lep == "elec";
  string name = isElec ? "el" : "mu";
  string outname = name + + "_metrel_sys";

  // Get the rate
  string plot = lep + "_qcdMC_all_metrel";
  TH1* temp = buildRate(m_mc.file, plot, m_mc.sname, "#slash{E}^{rel}_{T}", "", kBlack);

  // Get the normalizaiton
  string pnorm = lep + "_qcdMC_all_onebin";
  TH1* norm = buildRate(m_mc.file, pnorm, m_mc.sname, "onebin", "", kBlack);
  
  // Scale the rate
  temp->Scale(1/norm->GetBinContent(1));
  norm->Delete();

  TH1* rate = new TH1F(outname.c_str(),"#slash{E}^{rel}_{T} systematic",nMetbins, Metbins);
  rate->Sumw2();
  for(int bin = 1; bin<= rate->GetNbinsX(); ++bin){
    float bc = (temp->GetBinContent(bin) - 1);
    float be = temp->GetBinError(bin);
    rate->SetBinContent(bin, bc);
    rate->SetBinError(bin, be);
  }
  temp->Delete();
  return rate;
}
//------------------------------------------------------------//
TH1* FinalNewFake::getEtaSys(string lep)
{
  // Here we take the eta distribution and 
  // normalize this to the average fake rate 
  // (taken from one bin rate) and then take the differences from
  // 1 as the percentage.

  bool isElec = lep == "elec";
  string name = isElec ? "el" : "mu";
  string outname = name + + "_eta_sys";

  string plot = lep+"_qcdMC_all_l_eta_coarse";
  string pnorm = lep+"_qcdMC_all_onebin";
  
  // Get the rate and norm
  TH1* rate = buildRate(m_mc.file, plot, m_mc.sname, "|#eta|", "", kBlack);
  TH1* norm = buildRate(m_mc.file, pnorm, m_mc.sname, "onebin", "", kBlack);

  // Scale the rate
  rate->Scale(1/norm->GetBinContent(1));
  rate->SetName(outname.c_str());
  norm->Delete();
  for(int bin=1; bin<=rate->GetNbinsX(); ++bin){
    if( isElec ) rate->SetBinContent(bin, rate->GetBinContent(bin) - 1);
    else rate->SetBinContent(bin, 0); // mu consistent with 0
  }
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
//------------------------------------------------------------//
// Write to file
//------------------------------------------------------------//
void FinalNewFake::save(TH1* h)
{
  
  m_outfile->cd();
  h->Write();

}
//------------------------------------------------------------//
void FinalNewFake::save(TParameter<double> param)
{
  
  m_outfile->cd();
  param.Write();

}

//------------------------------------------------------------//
// Save the file
//------------------------------------------------------------//
void FinalNewFake::write()
{
  
  m_outfile->cd();
  m_outfile->Write();
  m_outfile->Close();

}

//------------------------------------------------------------//
// Simple plotting method
//------------------------------------------------------------//
void FinalNewFake::dumpPlot(TH1* hist)//, string save)
{
  TCanvas* c = makeCanvas("c");
  
  float xleg[] = {0.15, 0.3};
  float yleg[] = {0.85, 0.9};
  TLegend* leg = makeLegend(xleg,yleg);
  leg->SetTextSize(0.04);
  leg->AddEntry(hist);
  
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleOffset(0.85);
  hist->SetLabelSize(0.04, "Y");
  hist->SetLabelSize(0.04, "X");
  
  int bin = hist->GetMaximumBin();
  float max = hist->GetBinContent(bin) + hist->GetBinError(bin);
  hist->SetMaximum( max * 1.1 );

  hist->Draw("p");
  leg->Draw("same");
  
  string name = hist->GetName();
  string dir = basedir(m_outputfname);
  c->SaveAs((dir + name + ".eps").c_str());
  c->SaveAs((dir + name + ".png").c_str());
  delete c;
  
}
//------------------------------------------------------------//
// Dump stat error
//------------------------------------------------------------//
void FinalNewFake::dumpStat(TH1* hist)
{
  
  cout<<"Name: "<<hist->GetName();

  float stat = 0;
  for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
    float dev = (hist->GetBinError(bin)/hist->GetBinContent(bin)) * 100;
    if(stat < dev)
      stat = dev;
  }

  cout<<" Maximum deviation: "<<stat<<endl;

}

