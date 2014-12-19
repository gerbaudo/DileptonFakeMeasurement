// emacs -*- C++ -*-
#ifndef SUSY_WH_TUPLEMAKER_H
#define SUSY_WH_TUPLEMAKER_H

#include "SusyTest0/TupleMakerObjects.h"

#include <string>
#include <vector>

class TTree;
class TFile;
namespace Susy
{
class Lepton;
class Jet;
class Met;
}
// LeptonVector is defined in SusyDefs.h, but that's a huge include just for one def...refactor
typedef std::vector<Susy::Lepton*> LeptonVector;
typedef std::vector<Susy::Jet*>    JetVector;

namespace susy
{
namespace wh
{
/*!
  A class to save the information from SusyNt to a simpler ntuple.
  
  Details:
  This class is meant to create small ntupled for faster turnaround.
  The nutples store the information relative to the following objects:
  - two leading leptons
  - met
  - jets
  - other leptons
  - event variables
  The information is converted from SusyNt classes to smaller and
  simpler objects (see TupleMakerObjects.h)
  
  davide.gerbaudo@gmail.com
  November 2013
*/
class TupleMaker {
public:
    TupleMaker(const std::string &outFilename, const std::string &treename, bool delayInit=true);
    ~TupleMaker();
    bool init(const std::string &outFilename, const std::string &treename);
    bool close();
    bool fill(const double weight, const unsigned int run, const unsigned int event,
              const Susy::Lepton &l0, const Susy::Lepton &l1, const Susy::Met &met,
              const LeptonVector &otherLeptons, const JetVector &jets);
    //! methods to assign the pieces of info that are not accessible from Lepton (mostly fake-related)
    TupleMaker& setL0IsTight(bool v) { l0_.setIsTight(v); return *this; }
    TupleMaker& setL1IsTight(bool v) { l1_.setIsTight(v); return *this; }
    TupleMaker& setL0Source(int v) { l0_.setSource(v); return *this; }
    TupleMaker& setL1Source(int v) { l1_.setSource(v); return *this; }
    TupleMaker& setL0EtConeCorr(float v) { l0_.setEtConeCorr(v); return *this; }
    TupleMaker& setL1EtConeCorr(float v) { l1_.setEtConeCorr(v); return *this; }
    TupleMaker& setL0PtConeCorr(float v) { l0_.setPtConeCorr(v); return *this; }
    TupleMaker& setL1PtConeCorr(float v) { l1_.setPtConeCorr(v); return *this; }
    TupleMaker& setHas2ltrigmatch(const int &v) { eventPars_.setHas2ltrigmatch(v); return *this; }
    const TFile* file() const { return file_; }
    const TTree* tree() const { return tree_; }
    const FourMom& l0() const { return l0_; }
    const FourMom& l1() const { return l1_; }
    const FourMom& met() const { return met_; }
    const EventParameters& eventPars() const { return eventPars_; }
private: // rule of three 
    TupleMaker(const TupleMaker&);
    TupleMaker& operator=(const TupleMaker&);
private:
    bool initFile(const std::string &outFilename);
    bool initTree(const std::string &treename);
    bool initTreeBranches();
private:
    TFile *file_;
    TTree *tree_;
    FourMom l0_, l1_, met_;
    std::vector<FourMom> jets_, lowptLepts_;
    EventParameters eventPars_;
}; // end TupleMaker

} // namespace wh
} // namespace susy

#endif // end include guard
