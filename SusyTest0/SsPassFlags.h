#ifndef SsPassFlags_h
#define SsPassFlags_h

#include <string>

struct SsPassFlags {
  SsPassFlags() { reset(); }
  void reset() {
    eq2l = tauVeto = trig2l = trig2lmatch = true2l = sameSign = fjveto = bjveto = false;
    ge1j = eq1j = ge2j = false;
    lepPt = zllVeto = mtllmet = ht = metrel = veto3rdL = false;
  }
  bool passLpt() const {
    return (eq2l & tauVeto & trig2l & trig2lmatch & true2l & sameSign
            & fjveto & bjveto & ge1j & lepPt);
  }
  bool passAll() const { return (passLpt() & zllVeto & mtllmet & ht & metrel & veto3rdL); }
  bool eq2l, tauVeto, trig2l, trig2lmatch, true2l, sameSign, fjveto, bjveto;
  bool ge1j, eq1j, ge2j;
  bool lepPt, zllVeto, mtllmet, ht, metrel, veto3rdL;
  std::string str() const;
};

#endif
