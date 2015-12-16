// emacs -*- C++ -*-
#ifndef SUSY_WH_TUPLEMAKEROBJECTS_H
#define SUSY_WH_TUPLEMAKEROBJECTS_H

#include "SusyNtuple/SusyNt.h"

#include <iostream>
#include <vector>

namespace susy
{
namespace wh
{

using Susy::Lepton;
using Susy::Electron;
using Susy::Muon;
using Susy::Jet;
using Susy::Met;
struct FourMom {
    double px, py, pz, E;
    bool isMu, isEl, isJet;
    bool isTight;
    bool isTightPp; ///< used only for electron to compute tight/loose requirements offline
    int source; // see FakeLeptonSources
    unsigned int trigFlags;
    double charge, d0Signif, z0SinTheta, etCone, ptCone, mv1;
    double etConeCorr, ptConeCorr;
    FourMom() : px(0), py(0), pz(0), E(0),
                isMu(false), isEl(false), isJet(false),
                isTight(false),
                isTightPp(false),
                source(-1),
                trigFlags(0),
                charge(0), d0Signif(0), z0SinTheta(0), etCone(0), ptCone(0), mv1(0),
                etConeCorr(0), ptConeCorr(0) {}
#ifndef __CINT__
// cint is not able to parse 'complex' code; see
// http://root.cern.ch/drupal/content/interacting-shared-libraries-rootcint
    FourMom& set4mom(const Lepton &l) {
        const bool unbiased(true);
        px=l.Px(); py=l.Py(); pz=l.Pz(); E=l.E();
        charge = l.q;
        d0Signif = l.d0Sig(unbiased);
        z0SinTheta = l.z0SinTheta(unbiased);
        ptCone = l.ptcone30;
        return *this;
    }
    FourMom& set4mom(const Jet &j)    { px=j.Px(); py=j.Py(); pz=j.Pz(); E=j.E(); return *this; }
    FourMom& setMu(const Lepton &l) {
        isMu=true; isEl = isJet = false;
        if(const Muon* m = dynamic_cast<const Muon*>(&l)) etCone = m->etcone30;
        trigFlags = l.trigFlags;
        return set4mom(l);
    }
    FourMom& setEl(const Lepton &l) {
        isEl=true; isMu = isJet = false;
        if(const Electron *e = dynamic_cast<const Electron*>(&l)){
            etCone = e->topoEtcone30Corr;
            isTightPp = e->tightPP;
        }
        trigFlags = l.trigFlags;
        return set4mom(l);
    }
    FourMom& setIsTight(bool v) { isTight = v; return *this; }
    FourMom& setSource(int s) { source = s; return *this; }
    FourMom& setEtConeCorr(double v) { etConeCorr = v; return *this; }
    FourMom& setPtConeCorr(double v) { ptConeCorr = v; return *this; }
    FourMom& setJet(const Jet &j)   { isJet=true; isMu = isEl = false; mv1 = j.mv1; return set4mom(j); }
    FourMom& setMet(const Met &m)   { isJet=isMu=isEl=false; px=m.lv().Px(); py=m.lv().Py(); E=m.lv().E(); return *this; }
#endif // end ifndef CINT
}; // end FourMom

struct EventParameters {
    double weight;
    unsigned int eventNumber;
    unsigned int runNumber;
    int has2ltrigmatch;
    EventParameters() : weight(0), eventNumber(0), runNumber(0), has2ltrigmatch(0) {}
#ifndef __CINT__
    EventParameters& setWeight(const double &w) { weight=w; return *this; }
    EventParameters& setEvent(const unsigned int &e) { eventNumber=e; return *this; }
    EventParameters& setRun(const unsigned int &r) { runNumber=r; return *this; }
    EventParameters& setHas2ltrigmatch(const int &v) { has2ltrigmatch=v; return *this; }
#endif
};

} // namespace wh
} // namespace susy

#endif // end include guard
