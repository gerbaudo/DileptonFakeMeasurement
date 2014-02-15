#include "SusyTest0/SsPassFlags.h"
#include "SusyNtuple/SusyNt.h"
#include <sstream>      // std::ostringstream

//-----------------------------------------
SsPassFlags& SsPassFlags::updateLlFlags(const Susy::Lepton &l0, const Susy::Lepton &l1)
{
    bool e0(l0.isEle()), m0(l0.isMu());
    bool e1(l1.isEle()), m1(l1.isMu());
    ee = e0 && e1;
    em = ((e0 && m1) || (m0 && e1));
    mm = m0 && m1;
    return *this;
}
//-----------------------------------------
std::string SsPassFlags::str() const
{
  std::ostringstream oss;
  oss<<" eq2l: "<<eq2l
     <<" ll: "<<(ee ? "ee" : em ? "em" : mm ? "mm" : "??")
     <<" tauVeto: "<<tauVeto
     <<" trig2l: "<<trig2l
     <<" trig2lmatch: "<<trig2lmatch
     <<" true2l: "<<true2l
     <<" sameSign: "<<sameSign
     <<" veto3rdL: "<<veto3rdL
     <<" fjveto: "<<fjveto
     <<" bjveto: "<<bjveto
     <<" ge1j: "<<ge1j
     <<" eq1j: "<<eq1j
     <<" ge2j: "<<ge2j
     <<" lepPt: "<<lepPt
     <<" zllVeto: "<<zllVeto
     <<" dEtall: "<<dEtall
     <<" maxMt: "<<maxMt
     <<" mljj: "<<mljj
     <<" ht: "<<ht
     <<" metrel: "<<metrel
     <<" mtllmet: "<<mtllmet;
  return oss.str();
}
//-----------------------------------------
