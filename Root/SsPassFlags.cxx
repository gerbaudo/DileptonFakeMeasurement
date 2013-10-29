#include "SusyTest0/SsPassFlags.h"

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
     <<" fjveto: "<<fjveto
     <<" bjveto: "<<bjveto
     <<" ge1j: "<<ge1j
     <<" lepPt: "<<lepPt
     <<" zllVeto: "<<zllVeto
     <<" mtllmet: "<<mtllmet
     <<" ht: "<<ht
     <<" metrel: "<<metrel;
  return oss.str();
}
