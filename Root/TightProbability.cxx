#include "SusyTest0/TightProbability.h"

#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"

using std::cout;
using std::endl;

using Susy::TightProbability;


//----------------------------------------------------------
TightProbability::NumDenHisto::NumDenHisto(string name, int nbins, float min, float max):
  m_num((name+"_num").c_str(), (name+" numerator"  ).c_str(),nbins,min, max),
  m_den((name+"_den").c_str(), (name+" denominator").c_str(),nbins,min, max),
  m_min(0.0), m_max(0.0),
  m_widthFirst(0.0), m_widthLast(0.0)
{
  Sumw2();
  setMinMax();
}
TightProbability::NumDenHisto::NumDenHisto(string name, int nbins, const float* binEdges):
  m_num((name+"_num").c_str(), (name+" numerator"  ).c_str(),nbins, binEdges),
  m_den((name+"_den").c_str(), (name+" denominator").c_str(),nbins, binEdges),
  m_min(0.0), m_max(0.0),
  m_widthFirst(0.0), m_widthLast(0.0)
{
  Sumw2();
  setMinMax();
}
void TightProbability::NumDenHisto::Fill(bool alsoFillNum, float weight, float value)
{
  bool undeflow(value<m_min), overflow(value>m_max);
  value = (undeflow ? m_min+0.5*m_widthFirst : value);
  value = (overflow ? m_max-0.5*m_widthLast  : value);
  m_den.Fill(value, weight);
  if(alsoFillNum) m_num.Fill(value, weight);
}
void TightProbability::NumDenHisto::setMinMax() {
  if(const TAxis *ax = m_den.GetXaxis()){
    m_max = ax->GetXmax();
    m_min = ax->GetXmin();
    m_widthFirst = ax->GetBinWidth(1);
    m_widthLast = ax->GetBinWidth(ax->GetNbins());
  }
}
//----------------------------------------------------------
const float edgesFakePtbins[] = {5, 10, 15, 20, 25, 30, 35, 50, 100,500};
const int  nFakePtbins   = sizeof(edgesFakePtbins)/sizeof(float) - 1;
TightProbability::TightProbability():
  m_outFname("tightProbabilityOut.root"),
  m_outFile(0),
  m_isMC(false),
  m_evtWeight(0.0),
  m_h_pt_any  ("pt_any",   nFakePtbins, edgesFakePtbins),
  m_h_pt_real ("pt_real",  nFakePtbins, edgesFakePtbins),
  m_h_pt_hf   ("pt_hf",    nFakePtbins, edgesFakePtbins),
  m_h_pt_lf   ("pt_lf",    nFakePtbins, edgesFakePtbins),
  m_h_pt_conv ("pt_conv",  nFakePtbins, edgesFakePtbins),
  m_h_pt_mjet ("pt_mjet",  nFakePtbins, edgesFakePtbins),
  m_h_pt_other("pt_other", nFakePtbins, edgesFakePtbins),
  m_leptonOriginCounter(kLeptonOriginN, 0)
{}
//----------------------------------------------------------
TightProbability::~TightProbability() {}
//----------------------------------------------------------
void TightProbability::Begin(TTree *tree) {
  if(m_dbg>0) cout<<"MeasureFakeRate2::Begin"<<endl;
  SusySelection::Begin(0);
  initOutput(m_outFname);
}
//----------------------------------------------------------
void TightProbability::Terminate() {
  //if(m_dbg)
  dumpEventCounters();
  printLeptonOriginCounters();
  finalizeOutput();
}
//----------------------------------------------------------
Bool_t TightProbability::Process(Long64_t entry) {
  //  cout<<"entry "<<entry;
  GetEntry(entry);
  m_chainEntry++;
  clearObjects();
  m_ET = ET_Unknown;
  bool removeLepsFromIso(false);
  selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium);
  bool passEvent(selectEvent());
  //  cout<<(passEvent ? " passEvent " : "!passEvent\n");
  if(!passEvent) return true;
  m_ET = getDiLepEvtType(m_baseLeptons);
  bool passSr = passSrSs(m_ET, WH_SRSS1, m_baseLeptons, m_signalTaus, m_signalJets2Lep,
                         m_met);
  //  cout<<(passSr ? " passSrSs \n" : "!passSrSs\n");
  if(!passSr) return true;
  //  LeptonVector tags, probes;
  vector<float> pts(m_baseLeptons.size());
  for(size_t iL=0; iL<m_baseLeptons.size(); ++iL) {
    Lepton *l = m_baseLeptons[iL];
    bool fillNum = isSignalLepton(l, m_baseElectrons, m_baseMuons,
                                  nt.evt()->nVtx,nt.evt()->isMC);
    bool fn(fillNum);
    float weight(1.0), w(weight); // assume ~uniform weight; we're using a ratio...
    float pt(l->Pt());
    LeptonOrigin lo(getLeptonOrigin(l));
    m_h_pt_any.Fill(fillNum, weight, pt);
    vector<size_t> &oc = m_leptonOriginCounter;
    switch(lo){
    case kReal        : m_h_pt_real .Fill(fn, w, pt); oc[kReal         ]++; break;
    case kHeavyFlavor : m_h_pt_hf   .Fill(fn, w, pt); oc[kHeavyFlavor  ]++; break;
    case kLigthFlavor : m_h_pt_lf   .Fill(fn, w, pt); oc[kLigthFlavor  ]++; break;
    case kConversion  : m_h_pt_conv .Fill(fn, w, pt); oc[kConversion   ]++; break;
    case kMultijet    : m_h_pt_mjet .Fill(fn, w, pt); oc[kMultijet     ]++; break;
    default           : m_h_pt_other.Fill(fn, w, pt); oc[kUnknownOrigin]++; break;
    }
    pts[iL] = pt;
  } // end for(iL)
//   cout<<" pts: ";
//   copy(pts.begin(), pts.end(), ostream_iterator<float>(cout, " "));
//   cout<<endl;
  return true;
}
//----------------------------------------------------------
void TightProbability::initOutput(string outName) {
  if(m_dbg>0) cout<<"TightProbability::initOutput Creating file: "<<outName<<endl;
  m_outFile = new TFile(outName.c_str(),"recreate");
  m_outFile->cd();
  m_h_pt_any  .SetDirectory(m_outFile);
  m_h_pt_real .SetDirectory(m_outFile);
  m_h_pt_hf   .SetDirectory(m_outFile);
  m_h_pt_lf   .SetDirectory(m_outFile);
  m_h_pt_conv .SetDirectory(m_outFile);
  m_h_pt_mjet .SetDirectory(m_outFile);
  m_h_pt_other.SetDirectory(m_outFile);
}
//----------------------------------------------------------
void TightProbability::finalizeOutput() {
  if(m_dbg>0) cout<<"TightProbability::finalizeOutput"<<endl;
  if(m_dbg>1) m_outFile->ls();
  m_outFile->Write();
  m_outFile->Close();
  m_outFile->Delete();
}
//----------------------------------------------------------
TightProbability& TightProbability::setOutputFilename(const std::string &s) {
  m_outFname = s;
  return *this;
}
//----------------------------------------------------------
TightProbability::LeptonOrigin TightProbability::getLeptonOrigin(const Lepton* l){
  if( isRealLepton(l) ) return kReal;
  if( isHFLepton(l) )   return kHeavyFlavor;
  if( isLFLepton(l) )   return kLigthFlavor;
  if( isConvLepton(l) ) return kConversion;
  return kUnknownOrigin;
}
//----------------------------------------------------------
TightProbability::TightLoosePairType TightProbability::getTightLoosePairType(const Lepton* tag,
                                                           const Lepton* probe){
  if( !tag ) return kLooseLoose;
  // DG : this is not well defined for !probe.
  int nvtx  = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  bool ptSorted(tag->Pt() >= probe->Pt());
  const Lepton *l0 = ptSorted ? tag : probe;
  const Lepton *l1 = ptSorted ? probe : tag;
  bool l0isTight(isSignalLepton(l0, m_baseElectrons, m_baseMuons, nvtx, isMC));
  bool l1isTight(isSignalLepton(l1, m_baseElectrons, m_baseMuons, nvtx, isMC));
  return (l0isTight
          ? (l1isTight ? kTightTight : kTightLoose)
          : (l1isTight ? kLooseTight : kLooseLoose));
}
//----------------------------------------------------------
float TightProbability::getPtcone(const Lepton* lep){
  if( lep->isEle() )  return lep->ptcone30;
  else  return ((Muon*) lep)->ptcone30ElStyle;
}
//----------------------------------------------------------
float TightProbability::getEtcone(const Lepton* lep){
  if( lep->isEle() ) return ((Electron*) lep)->topoEtcone30Corr;
  return ((Muon*) lep)->etcone30;
}
//----------------------------------------------------------
void TightProbability::printLeptonOriginCounters() const {
    const vector<size_t> &oc = m_leptonOriginCounter;
    cout<<"TightProbability: leptonOriginCounters"<<endl
        <<"kReal          : "<<oc[kReal         ]<<endl
        <<"kHeavyFlavor   : "<<oc[kHeavyFlavor  ]<<endl
        <<"kLigthFlavor   : "<<oc[kLigthFlavor  ]<<endl
        <<"kConversion    : "<<oc[kConversion   ]<<endl
        <<"kMultijet      : "<<oc[kMultijet     ]<<endl
        <<"kUnknownOrigin : "<<oc[kUnknownOrigin]<<endl
        <<endl;
}
//----------------------------------------------------------
