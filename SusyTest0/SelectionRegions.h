#ifndef SUSYWH_SELECTIONREGIONS_H
#define SUSYWH_SELECTIONREGIONS_H

// enumFactory provides the functions GetString and SelectionRegions2str
#include "SusyTest0/enumFactory.h"

namespace susywh
{

#define SOME_ENUM(VALS) \
  VALS( kNoSel         ,) \
  VALS( kSr6base       ,) \
  VALS( kSr6           ,) \
  VALS( kSr7base       ,) \
  VALS( kSr7Nj         ,) \
  VALS( kSr7NjZttVeto  ,) \
  VALS( kSr7NjPtTot    ,) \
  VALS( kSr7NjMll      ,) \
  VALS( kSr7           ,) \
  VALS( kSr8base       ,) \
  VALS( kCr8lpt        ,) \
  VALS( kCr8ee         ,) \
  VALS( kCr8mm         ,) \
  VALS( kSr8           ,) \
  VALS( kSr9base       ,) \
  VALS( kCr9lpt        ,) \
  VALS( kSr9           ,) \
  VALS( kN             ,) \

  DECLARE_ENUM(SelectionRegions, SOME_ENUM)

} //end namespace susywh

#endif // end include guard
