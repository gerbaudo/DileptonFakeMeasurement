#include "SusyTest0/kinematic.h"

#include "SusyNtuple/SusyNtTools.h"

#include "TLorentzVector.h"

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
        v.isEm = (isEl0 && isMu0) || (isMu0 && isEl1);
        v.isMm = isMu0 && isMu1;
        v.pt0 = l0.Pt();
        v.pt1 = l1.Pt();
        v.mll = (l0+l1).M();
        LeptonVector lepts;
        lepts.push_back(&l0);
        lepts.push_back(&l1);
        v.metrel = SusyNtTools::getMetRel(met, lepts, jets);
        if     (v.numCentralLightJets==1) v.mlj  = swk::mlj (l0, l1, *jets[0]);
        else if(v.numCentralLightJets >1) v.mljj = swk::mljj(l0, l1, *jets[0], *jets[1]);
        v.mt0 = 0.0; //transverseMass(ll, met->lv()) // move from criteria to kin
        v.mt1 = 0.0; //transverseMass(ll, met->lv()) // move from criteria to kin
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

