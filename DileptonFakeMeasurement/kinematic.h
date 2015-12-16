#ifndef SUSY_WH_KIN_H
#define SUSY_WH_KIN_H

/*
  Functions to compute kinematic variables for WH events.

  davide.gerbaudo@gmail.com
  Jan 2014
*/

#include <cstddef> // size_t

#include "SusyNtuple/SusyNt.h" // Lepton, Jet, Met, and all that
#include "SusyNtuple/SusyDefs.h" // LeptonVector, JetVector and all that

namespace susy
{
namespace wh
{
namespace kin
{
  struct DilepVars { //! container for dilepton variables we want to compute once per event
    DilepVars() { reset(); }
    bool isEe, isEm, isMm;
    size_t numCentralLightJets;
    float pt0, pt1;
    float mll, detall;
    float metrel;
    float mlj, mljj;
    float mt0, mt1, mtllmet;
    float ht;
    bool l3veto;
    float mtmin() const { return mt0<mt1 ? mt0 : mt1; }
    float mtmax() const { return mt0>mt1 ? mt0 : mt1; }
    void reset() {
      isEe = isEm = isMm = false;
      numCentralLightJets = 0;
      pt0 = pt1 =  mll = detall = metrel = mlj = mljj = mt0 = mt1 = mtllmet = 0.0;
      ht = 0.0;
      l3veto = false;
    }
  };
  DilepVars compute2lVars(const LeptonVector &leptons, const Susy::Met *met, const JetVector &jets);
  float mlj (const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &j);
  float mljj(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &j0,const TLorentzVector &j1);
  float transverseMass(const TLorentzVector &lep, const TLorentzVector &met);
  //! DG : static copy of SusyNtTools::Meff, which uses all leptons, all jets, and met; is this what we want?
  float meff(const TLorentzVector &l0, const TLorentzVector &l1, const Susy::Met* met, const JetVector &jets);

} // end kin
} // end wh
} // end susy

#endif
