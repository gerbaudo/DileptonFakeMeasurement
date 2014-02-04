#include "SusyTest0/kinematic.h"

#include "SusyNtuple/SusyNtTools.h"

#include "TLorentzVector.h"

#include <cmath> // sqrt
#include <math.h> // cos, fabs
#include <cassert>

namespace swk = susy::wh::kin;

swk::DilepVars swk::compute2lVars(const LeptonVector &leptons, const Susy::Met *met, const JetVector &jets)
{
    DilepVars v;
    v.numCentralLightJets = jets.size();
    if(leptons.size()>1) {
        Susy::Lepton &l0 = *leptons[0];
        Susy::Lepton &l1 = *leptons[1];
        bool isEl0(l0.isEle()), isEl1(l1.isEle()), isMu0(l0.isMu()), isMu1(l1.isMu());
        assert(isEl0!=isMu0 && isEl1!=isMu1); // assuming we're only dealing with electrons or muons
        v.isEe = isEl0 && isEl1;
        v.isEm = ((isEl0 && isMu1) ||
                  (isMu0 && isEl1)  );
        v.isMm = isMu0 && isMu1;
        v.pt0 = l0.Pt();
        v.pt1 = l1.Pt();
        TLorentzVector ll(l0+l1);
        v.mll = ll.M();
        v.detall = fabs(l0.Eta() - l1.Eta());
        LeptonVector lepts;
        lepts.push_back(&l0);
        lepts.push_back(&l1);
        v.metrel = SusyNtTools::getMetRel(met, lepts, jets);
        if     (v.numCentralLightJets==1) v.mlj  = swk::mlj (l0, l1, *jets[0]);
        else if(v.numCentralLightJets >1) v.mljj = swk::mljj(l0, l1, *jets[0], *jets[1]);
        v.mt0 = swk::transverseMass(l0, met->lv());
        v.mt1 = swk::transverseMass(l1, met->lv());
        v.ht = swk::meff(l0, l1, met, jets);
        v.mtllmet = transverseMass(ll, met->lv());
        // NB l3veto calculation should be moved here -- DG 2014-01-29
    }
    return v;
}

float swk::mlj(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &j)
{
    float dr0(j.DeltaR(l0)), dr1(j.DeltaR(l1));
    return dr0<dr1 ? (j+l0).M() : (j+l1).M();
}

float swk::mljj(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &j0, const TLorentzVector &j1)
{
    return swk::mlj(l0, l1, j0+j1);
}
//-----------------------------------------
float swk::transverseMass(const TLorentzVector &lep, const TLorentzVector &met)
{
  return std::sqrt(2.0 * lep.Pt() * met.Et() *(1-cos(lep.DeltaPhi(met))) );
}
//-----------------------------------------
float swk::meff(const TLorentzVector &l0, const TLorentzVector &l1, const Susy::Met* met, const JetVector &jets)
{
  float meff = 0;
  meff += l0.Pt();
  meff += l1.Pt();
  for(size_t i=0; i<jets.size(); i++){ if(jets[i]->Pt() > 20.0) meff += jets[i]->Pt(); }
  meff += met->Et;
  return meff;
}
//-----------------------------------------
