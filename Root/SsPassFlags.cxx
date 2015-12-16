#include "DileptonFakeMeasurement/SsPassFlags.h"

#include <sstream>      // std::ostringstream

//-----------------------------------------
std::string SsPassFlags::str() const
{
  std::ostringstream oss;
  oss<<" eq2l: "<<eq2l
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
