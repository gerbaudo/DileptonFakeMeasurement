#ifndef SUSY_CRITERIA_H
#define SUSY_CRITERIA_H

/*
  Functions defining criteria to select objects and events.

  davide.gerbaudo@gmail.com
  Aug 2013
 */

#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyDefs.h"

namespace susy
{

bool isRealLepton(const Susy::Lepton* lep);
bool isFakeLepton(const Susy::Lepton* lep);
bool isConvLepton(const Susy::Lepton* lep);
bool isHFLepton(const Susy::Lepton* lep);
bool isLFLepton(const Susy::Lepton* lep);
bool isTrueDilepton(const LeptonVector &leptons);
bool passEleD0S(const LeptonVector &leptons, float maxVal);

} // end namespace susy

#endif
