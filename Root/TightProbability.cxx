#include "SusyTest0/TightProbability.h"

#include <iostream>

#include "TFile.h"
#include "TH1F.h"

using std::cout;
using std::endl;

using Susy::TightProbability;

//----------------------------------------------------------
TightProbability::TightProbability():
  m_outFname("tightProbabilityOut.root"),
  m_outFile(0),
  m_isMC(false),
  m_evtWeight(0.0),
  m_AltIso(false),
  m_ch(ET_Unknown)
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
  finalizeOutput();
}
//----------------------------------------------------------
Bool_t TightProbability::Process(Long64_t entry) {return true;}
//----------------------------------------------------------
void TightProbability::initOutput(string outName) {
  if(m_dbg>0) cout<<"TightProbability::initOutput Creating file: "<<outName<<endl;
  m_outFile = new TFile(outName.c_str(),"recreate");
  m_outFile->cd();
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
