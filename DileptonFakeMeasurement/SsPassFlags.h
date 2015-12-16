#ifndef SsPassFlags_h
#define SsPassFlags_h

#include <string>

struct SsPassFlags {
  SsPassFlags() { reset(); }
  void reset() {
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
  bool eq2l, tauVeto, trig2l, trig2lmatch, true2l, sameSign, veto3rdL, fjveto, bjveto;
  bool ge1j, eq1j, ge2j;
  bool lepPt, zllVeto, dEtall, maxMt, mljj, ht, metrel, mtllmet;
  std::string str() const;
};

#endif
