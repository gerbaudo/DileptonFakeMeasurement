#ifndef SUSY_CRITERIA_H
#define SUSY_CRITERIA_H

/*
  Functions defining criteria to select objects and events.

  \todo: move it to a namespace

  davide.gerbaudo@gmail.com
  Aug 2013
 */

#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyDefs.h"

//#include "SusyNtuple/SusyNtAna.h"
//#include "SusyNtuple/DilTrigLogic.h"
//#include "SusyNtuple/SusyNtTools.h"
//#include "SusyNtuple/SusyDefs.h"

//#include "SUSYTools/SUSYObjDef.h"



    // Idendification methods
bool isRealLepton(const Susy::Lepton* lep);
bool isFakeLepton(const Susy::Lepton* lep);
bool isConvLepton(const Susy::Lepton* lep);
bool isHFLepton(const Susy::Lepton* lep);
bool isLFLepton(const Susy::Lepton* lep);
bool isTrueDilepton(const LeptonVector &leptons);
bool passEleD0S(const LeptonVector &leptons, float maxVal);


#endif
