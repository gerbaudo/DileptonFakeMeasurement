#include "SusyTest0/TightProbability.h"

using Susy::TightProbability;

TightProbability::TightProbability() {}
TightProbability::~TightProbability() {}
void TightProbability::Begin(TTree *tree) {}
void TightProbability::Terminate() {}
Bool_t TightProbability::Process(Long64_t entry) {return true;}
void TightProbability::initHistos(string outName) {}

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
