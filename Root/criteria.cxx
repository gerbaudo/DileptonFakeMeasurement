#include "SusyTest0/criteria.h"

#include <cassert>
#include <math.h>   // cos
#include <numeric>  // std::accumulate

#include "Mt2/mt2_bisect.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/WhTruthExtractor.h"
#include "LeptonTruthTools/RecoTruthMatch.h"
#include "SusyTest0/DileptonAnalyticalSolver.h"

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
//----------------------------------------------------------
bool passZllVeto(const LeptonVector& l, float mllLo, float mllHi)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  float mll((*l[0] + *l[1]).M());
  return (mll<mllLo || mllHi<mll);
}
//----------------------------------------------------------
void swap(float &a, float &b) { float c(a); a=b; b=c; };
bool pass2LepPt(const LeptonVector& l, float minPt0, float minPt1)
{
  if(l.size()<2) return false;
  float pt0(l[0]->Pt()), pt1(l[1]->Pt());
  if(pt0 < pt1) swap(pt0, pt1);
  return (pt0>minPt0 && pt1>minPt1);
}
//----------------------------------------------------------
bool passPtllMin(const LeptonVector& l, float minPt)
{
  if(l.size()<2) return false;
  return TLorentzVector(*l[0] + *l[1]).Pt() > minPt;
}
//----------------------------------------------------------
bool passPtTot(const LeptonVector& l, const JetVector& j, const Susy::Met* m, float maxPtTot)
{
  if(!m) return false;
  TLorentzVector mlv(m->lv()), ll, jj;
  if (l.size()>0 && l[0]) ll += *l[0];
  if (l.size()>1 && l[1]) ll += *l[1];
  if (j.size()>0 && j[0]) jj += *j[0];
  if (j.size()>1 && j[1]) jj += *j[1];
  return (ll+jj+mlv).Pt() < maxPtTot;
}
//----------------------------------------------------------
bool passMllMax(const LeptonVector& l, float maxMll)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  return TLorentzVector(*l[0] + *l[1]).M() < maxMll;
}
//----------------------------------------------------------
bool passMllMin(const LeptonVector& l, float minVal)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  return minVal < TLorentzVector(*l[0] + *l[1]).M();
}
//----------------------------------------------------------
bool passDrllMax(const LeptonVector& l, float maxDr)
{
  if(l.size()<2) return false;
  return  l[0]->DeltaR(*(static_cast<TLorentzVector*>(l[1]))) < maxDr;
}


//-----------------------------------------
bool passdPhi(TLorentzVector v0, TLorentzVector v1, float cut)
{
  return v0.DeltaPhi(v1) > cut;
}
//-----------------------------------------
bool passMtLlMetMin(const LeptonVector& l, const Susy::Met* met, float minVal)
{
  if(l.size() < 2 || !l[0] || !l[1]) return false;
  TLorentzVector ll = (*l[0] + *l[1]);
  float mww=transverseMass(ll, met->lv());
  return (minVal < mww);
}
//----------------------------------------------------------
bool passMtMinlmetMin(const LeptonVector& l, const Susy::Met* met, float minVal)
{
  if( l.size() < 2 ) return false;
  const TLorentzVector *l0 = static_cast<TLorentzVector*>(l[0]);
  const TLorentzVector *l1 = static_cast<TLorentzVector*>(l[1]);
  float mtL0Met = transverseMass(*l0, met->lv());
  float mtL1Met = transverseMass(*l1, met->lv());
  return (mtL0Met < mtL1Met ? mtL0Met : mtL1Met) > minVal;
}
//----------------------------------------------------------
float computeMt2(const TLorentzVector &l0, const TLorentzVector &l1, const TLorentzVector &met)
{
  double pTMiss[3] = {0.0, met.Px(), met.Py()};
  double pA[3]     = {0.0, l0.Px(), l0.Py()};
  double pB[3]     = {0.0, l1.Px(), l1.Py()};
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pA,pB,pTMiss);
  mt2_event.set_mn(0); // LSP mass = 0 is Generic
  return mt2_event.get_mt2();
}
//----------------------------------------------------------
bool passMT2(const LeptonVector& leptons, const Susy::Met* met, float cut)
{
  if( leptons.size() < 2 ) return false;
  TLorentzVector metlv = met->lv();
  TLorentzVector l0    = *leptons.at(0);
  TLorentzVector l1    = *leptons.at(1);
  return (computeMt2(l0, l1, metlv) > cut);
}
//----------------------------------------------------------
bool passHtMin(const LeptonVector& leptons,
               const JetVector &jets,
               const Susy::Met* met,
               float minVal)
{ // DG : static copy of SusyNtTools::Meff, which uses all leptons, all jets, and met; is this what we want?
  float meff = 0;
  for(uint i=0; i<leptons.size(); i++) meff += leptons[i]->Pt();
  for(uint i=0; i<jets.size(); i++){
    if(jets[i]->Pt() > 20.0) meff += jets[i]->Pt();
  }
  meff += met->Et;
  return (minVal < meff);
}
//----------------------------------------------------------
bool passNlepMin(const LeptonVector &leptons, size_t minVal)
{
  // return (leptons.size() >= minVal);
  // DG we should define a m_signalLeptons2L instead
  // similar to Anyes' implementation
  size_t nLep=0;
  for(size_t i=0;i<leptons.size(); ++i){
    if(const Susy::Lepton* l = leptons[i])
      if(l->isMu() && fabs(l->Eta())>2.4) // 2L muon trig, see Anyes' email 2013-08-02
        return false;
    nLep++;
  }
  return nLep>=minVal;
}
//----------------------------------------------------------
bool passZtautauVeto(const LeptonVector& l, const JetVector& j, const Susy::Met* m, float widthZpeak)
{
  if(l.size()<2 || !l[0] || !l[1]) return false;
  if(!m)         return false;
  float mZ0(91.);
  return abs(mZTauTau(*l[0], *l[1], m->lv()) - mZ0) > widthZpeak;
}
//-----------------------------------------
float getLeptonEff2Lep(const LeptonVector &leptons)
{
  assert(leptons.size()>1);
  return leptons[0]->effSF * leptons[1]->effSF;
}
//-----------------------------------------
int pdgIdFromLep(const Lepton *l)
{
  // particles have positive codes, see doi:10.1146/annurev.ns.25.120175.003011
  int kPel(+11), kAel(-11), kPmu(+13), kAmu(-13), kUnknown(0);
  if     (l->isEle()) return (l->q < 0 ? kPel : kAel);
  else if(l->isMu() ) return (l->q < 0 ? kPmu : kAmu);
  else                return kUnknown;
}
//-----------------------------------------
float transverseMass(const TLorentzVector &lep, const TLorentzVector &met)
{
  return std::sqrt(2.0 * lep.Pt() * met.Et() *(1-cos(lep.DeltaPhi(met))) );
}
//-----------------------------------------
float mtWW(const TLorentzVector &ll, const TLorentzVector &met)
{
  using std::sqrt;
  float dphi = acos(cos(ll.Phi() - met.Phi()));
  float mll(ll.M()), ptll(ll.Pt()), ptvv(met.Pt());
  float mvv = (/*mvvTrue*/ true ? mll : 0.0);
  return sqrt(mll*mll + mvv*mvv + 2.0*(sqrt  (ptll*ptll + mll*mll)
                                       * sqrt(ptvv*ptvv + mvv*mvv)
                                       - ptll * ptvv * cos(dphi)));
}
//-----------------------------------------
float sumCosDeltaPhi(const TLorentzVector &l0, const TLorentzVector &l1,
                     const TLorentzVector &met)
{
  return cos(l0.Phi() - met.Phi()) + cos(l1.Phi() - met.Phi());
}
//-----------------------------------------
float addJetPt(float totPt, const Susy::Jet *j) { return totPt + j->Pt(); }
float sumEtEtMiss(const TLorentzVector &el, const TLorentzVector &mu,
			       const JetVector &jets, const TLorentzVector &met)
{
  return
    el.Et()
    + mu.Pt()
    + met.Et()
    + std::accumulate(jets.begin(), jets.end(), float(0.0), addJetPt);
}
//-----------------------------------------
int numberOfNeutrinoSolutions(const TLorentzVector &lPos, const TLorentzVector &lNeg,
                              const Susy::Jet &jet0, const Susy::Jet &jet1,
                              const TLorentzVector &met)
{
    double mWp=80.41, mWm=80.41;
    double mnu=0.0,   mnub=0.0;
    double mt=172.9,  mtb=172.9;
    double ETmiss[2] = {met.Px(), met.Py()};
    double b[4]  = {jet0.E(), jet0.Px(), jet0.Py(), jet0.Pz()};
    double bb[4] = {jet1.E(), jet1.Px(), jet1.Py(), jet1.Pz()};
    double lp[4] = {lPos.E(), lPos.Px(), lPos.Py(), lPos.Pz()};
    double lm[4] = {lNeg.E(), lNeg.Px(), lNeg.Py(), lNeg.Pz()};
    std::vector<double> pnux, pnuy, pnuz, pnubx, pnuby, pnubz, cd_diff;
    int cubic_single_root_cmplx;
    llsolver::DileptonAnalyticalSolver slv;
    slv.solve(ETmiss, b, bb, lp, lm, mWp, mWm, mt, mtb, mnu, mnub,
	      &pnux, &pnuy, &pnuz, &pnubx, &pnuby, &pnubz,
	      &cd_diff, cubic_single_root_cmplx);
    return pnubx.size();
}
//-----------------------------------------
float mZTauTau(const TLorentzVector &l0, const TLorentzVector &l1,
               const TLorentzVector &met)
{ // re-written based on HWWlvlvCode::calculate_METBasedVariables
  float px0(l0.Px()), py0(l0.Py());
  float px1(l1.Px()), py1(l1.Py());
  float pxm(met.Px()), pym(met.Py());
  float num( px0*py1 - py0*px1 );
  float den1( py1*pxm - px1*pym + px0*py1 - py0*px1 );
  float den2( px0*pym - py0*pxm + px0*py1 - py0*px1 );
  float x1 = ( den1 != 0.0  ? (num/den1) : 0.0);
  float x2 = ( den2 != 0.0  ? (num/den2) : 0.0);
  bool kinematicallyPossible(x1*x2 > 0.0);
  return (kinematicallyPossible ? (l0+l1).M() / std::sqrt(x1*x2) : -1.0);
}


} // end namespace susy
