#ifndef SUSY_CRITERIA_H
#define SUSY_CRITERIA_H

/*
  Functions defining criteria to select objects and events.

  davide.gerbaudo@gmail.com
  Aug 2013
 */

#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyDefs.h"

class TLorentzVector;

namespace susy
{
  bool isRealLepton(const Susy::Lepton* lep);
  bool isFakeLepton(const Susy::Lepton* lep);
  bool isConvLepton(const Susy::Lepton* lep);
  bool isHFLepton(const Susy::Lepton* lep);
  bool isLFLepton(const Susy::Lepton* lep);
  bool isTrueDilepton(const LeptonVector &leptons);
  bool passEleD0S(const LeptonVector &leptons, float maxVal);
  bool sameFlavor(const LeptonVector& leptons);
  bool oppositeFlavor(const LeptonVector& leptons);
  bool sameSign(const LeptonVector& leptons);
  bool oppositeSign(const LeptonVector& leptons);
  bool passHtautauVeto(int hdecay);
  bool passZllVeto(const LeptonVector& l, float mllLo, float mllHi);
  bool pass2LepPt(const LeptonVector& l, float minPt0, float minPt1);
  bool passPtllMin(const LeptonVector& l, float minPt);
  bool passPtTot(const LeptonVector& l, const JetVector& j, const Susy::Met* m, float maxPtTot);
  bool passMllMax(const LeptonVector& l, float maxMll);
  bool passMllMin(const LeptonVector& l, float minVal);
  bool passDrllMax(const LeptonVector& l, float maxDr);
  bool passdPhi(TLorentzVector v0, TLorentzVector v1, float cut);
  bool passMtLlMetMin(const LeptonVector& leptons, const Susy::Met* met, float minVal=50.0);
  bool passMtMinlmetMin(const LeptonVector& leptons, const Susy::Met* met, float minVal=50.0);
  float computeMt2(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &met);
  bool passMT2(const LeptonVector& leptons, const Susy::Met* met, float cut);
  bool passHtMin(const LeptonVector& l, const JetVector &j, const Susy::Met* met, float minVal);
  bool passNlepMin(const LeptonVector &leptons, size_t minVal);
  bool passZtautauVeto(const LeptonVector& l, const JetVector& j, const Susy::Met* m, float widthZpeak=40.0);
  float getLeptonEff2Lep(const LeptonVector &leptons);
  int pdgIdFromLep(const Susy::Lepton *l);
  float transverseMass(const TLorentzVector &lep, const TLorentzVector &met);
  //! redundant? mt with non-zero m_vv? see https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/ataffard/SusyWeakProdAna/trunk/Root/PhysicsTools.cxx
  float mtWW(const TLorentzVector &ll, const TLorentzVector &met);
  //! compute tau-tau mass assuming that v's are collinear with leptons and responsible for all MET
  float mZTauTau(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &met);
  //! \f$ \Sum cos \Delta\phi \f$ used in CERN-PH-EP-2011-097
  float sumCosDeltaPhi(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &met);
  //! \f$ \Sum E_{T} + E_{T}^{miss} \f$ used in CERN-PH-EP-2011-097
  float sumEtEtMiss(const TLorentzVector &el, const TLorentzVector &mu,
                    const JetVector &jets, const TLorentzVector &met);
  //! compute the number of kinematically allowed neutrino solutions
  /*! This code is based on hep-ph/0603011. The number of possible
    solutions ranges from 0 to 4. However, you should consider
    calling the function twice (swapping the two jets), unless you
    can tell which one is the b and which one is the bbar.
  */
  int numberOfNeutrinoSolutions(const TLorentzVector &lPos, const TLorentzVector &lNeg,
                                const Susy::Jet &jet0, const Susy::Jet &jet1,
                                const TLorentzVector &met);
} // end namespace susy

#endif
