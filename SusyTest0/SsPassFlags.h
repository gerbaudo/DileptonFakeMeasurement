#ifndef SsPassFlags_h
#define SsPassFlags_h

#include <string>

namespace Susy
{
class Lepton;
}

struct SsPassFlags {
  SsPassFlags() { reset(); }
  void reset() {
    ee = em = mm = false;
    eq2l = tauVeto = trig2l = trig2lmatch = true2l = sameSign = veto3rdL = fjveto = bjveto = false;
    ge1j = eq1j = ge2j = false;
    lepPt = zllVeto = dEtall = maxMt = mljj = ht = metrel = mtllmet = false;
  }
  bool passLpt() const {
    return (eq2l & tauVeto & trig2l & trig2lmatch & true2l & sameSign
            & veto3rdL & fjveto & bjveto & ge1j & lepPt);
  }
  bool passCommonCriteria() const {
    return (eq2l & tauVeto & trig2l & trig2lmatch & true2l & sameSign
            & veto3rdL & fjveto & bjveto & ge1j);
  }
  bool passKinCriteria() const { //!< all the kin criteria that depend on (ee/em/mm) x (1j/23j); see SusySelection::passSrWh*j()
    return veto3rdL && lepPt && zllVeto && dEtall && maxMt && mljj && ht && metrel && mtllmet;
  }
  SsPassFlags& updateLlFlags(const Susy::Lepton &l0, const Susy::Lepton &l1);
  bool ee, em, mm;
  bool eq2l, tauVeto, trig2l, trig2lmatch, true2l, sameSign, veto3rdL, fjveto, bjveto;
  bool ge1j, eq1j, ge2j;
  bool lepPt, zllVeto, dEtall, maxMt, mljj, ht, metrel, mtllmet;
  std::string str() const;
};

#endif
