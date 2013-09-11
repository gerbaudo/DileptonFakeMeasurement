#include "SusyTest0/criteria.h"


#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/WhTruthExtractor.h"
#include "LeptonTruthTools/RecoTruthMatch.h"

using Susy::Lepton;

namespace susy
{
//-----------------------------------------
bool isRealLepton(const Lepton* lep)
{
  // Updated way of handling real and fake leptons using LeptonTruthTools
  return (lep->truthType == RecoTruthMatch::PROMPT);
}
//-----------------------------------------
bool isFakeLepton(const Lepton* lep)
{
  return !isRealLepton(lep);
}
//-----------------------------------------
bool isConvLepton(const Lepton* lep)
{
  bool isConverted(lep->truthType == RecoTruthMatch::CONV);
  bool isChargeFlip(lep->isEle() ?
                    static_cast<const Susy::Electron*>(lep)->isChargeFlip : false);
  return isConverted && !isChargeFlip;
}
//-----------------------------------------
bool isHFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::HF);
}
//-----------------------------------------
bool isLFLepton(const Lepton* lep)
{
  return (lep->truthType == RecoTruthMatch::LF);
}
//-----------------------------------------
bool isTrueDilepton(const LeptonVector &leptons)
{
  // Maybe not 100% kosher, but I just want to make sure I am not
  // double counting, so I require dilepton events to be real
  if( leptons.size() != 2 ) return false;
  bool l0_real = isRealLepton(leptons[0]);
  bool l1_real = isRealLepton(leptons[1]);
  bool l0_cf = leptons[0]->isEle() ? ((Susy::Electron*) leptons[0])->isChargeFlip : false;
  bool l1_cf = leptons[1]->isEle() ? ((Susy::Electron*) leptons[1])->isChargeFlip : false;
  return l0_real && l1_real && !l0_cf && !l1_cf; // ignoring charge flip
}
//-----------------------------------------
bool passEleD0S(const LeptonVector &leptons, float maxVal)
{
  for(size_t i=0; i<leptons.size(); ++i){
      const Susy::Lepton* l = leptons[i];
      if(l->isEle() && (fabs(l->d0Sig(true)) > maxVal)) return false;
  } // end for(i)
  return true;
}
//-----------------------------------------
bool sameFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return (leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
//-----------------------------------------
bool oppositeFlavor(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return !(leptons.at(0)->isMu() == leptons.at(1)->isMu());
}
//-----------------------------------------
bool sameSign(const LeptonVector& leptons)
{
  if( leptons.size() < 2 ) return false;
  return leptons.at(0)->q * leptons.at(1)->q > 0;
}
//-----------------------------------------
bool oppositeSign(const LeptonVector& leptons)
{
  return !(sameSign(leptons));
}
//-----------------------------------------
bool passHtautauVeto(int hdecay)
{
  return (hdecay!=WhTruthExtractor::kPtauAtau);
}


} // end namespace susy
