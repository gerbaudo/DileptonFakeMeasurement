#ifndef SUSYWH_SELECTIONREGIONS_H
#define SUSYWH_SELECTIONREGIONS_H

// enumFactory provides the functions GetString and SelectionRegions2str
#include "DileptonFakeMeasurement/enumFactory.h"

namespace susywh
{

#define SOME_ENUM(VALS) \
  VALS( kNoSel           ,) \
  VALS( kSr6base         ,) \
  VALS( kSr6             ,) \
  VALS( kSr7base         ,) \
  VALS( kSr7Nj           ,) \
  VALS( kSr7NjZttVeto    ,) \
  VALS( kSr7NjPtTot      ,) \
  VALS( kSr7NjMll        ,) \
  VALS( kSr7             ,) \
  VALS( kSr8base         ,) \
  VALS( kCr8lpt          ,) \
  VALS( kCr8ee           ,) \
  VALS( kCr8mm           ,) \
  VALS( kSr8             ,) \
  VALS( kSr9base         ,) \
  VALS( kCr9lpt          ,) \
  VALS( kSr9             ,) \
  VALS( kN               ,) \
  /* DG These are needed to compute the SF and the rates (pseudo t&p) */ \
  VALS( kCR_Real         ,)  /* Real Z window */           \
  VALS( kCR_SideLow      ,)  /* Real Side band Lower */    \
  VALS( kCR_SideHigh     ,)  /* Real Side band higher */   \
  VALS( kCR_HF           ,)  /* HF b-bbar tag and probe */ \
  VALS( kCR_HF_high      ,)  /* HF b-bbar t&p w/ met<80 */ \
  VALS( kCR_LFZjet       ,)  /* LF Z+jet tag and probe */  \
  VALS( kCR_LFWjet       ,)  /* LF W+jet tag and probe */  \
  VALS( kCR_Conv         ,)  /* Zmumu tag and probe */     \
  VALS( kCR_CFlip        ,)  /* Charge Flip */             \
  VALS( kCR_MCHeavy      ,)  /* MC cr heavy flavor */      \
  VALS( kCR_MCLight      ,)  /* MC cr light flavor */      \
  VALS( kCR_MCConv       ,)  /* MC cr conversions */       \
  VALS( kCR_MCQCD        ,)  /* MC cr QCD = LF+HF */       \
  VALS( kCR_MCALL        ,)  /* MC All */                  \
  VALS( kCR_MCReal       ,)  /* MC cr Real */              \
  VALS( kCR_MCNone       ,)  /* MC cr w/out metrel cuts */ \
  /* DG These are the SR, needed to compute the compositions. */ \
  VALS( kCR_SRmT2a       ,) \
  VALS( kCR_SRmT2b       ,) \
  VALS( kCR_SRmT2c       ,) \
  VALS( kCR_SRWWa        ,) \
  VALS( kCR_SRWWb        ,) \
  VALS( kCR_SRWWc        ,) \
  VALS( kCR_SRZjets      ,) \
  VALS( kCR_VRSS         ,) \
  VALS( kCR_CRWWMet      ,) \
  VALS( kCR_CRWWmT2      ,) \
  VALS( kCR_CRTopMet     ,) \
  VALS( kCR_CRTopmT2     ,) \
  VALS( kCR_CRZVMet      ,) \
  VALS( kCR_CRZVmT2_90   ,) \
  VALS( kCR_CRZVmT2_120  ,) \
  VALS( kCR_CRZVmT2_150  ,) \
  VALS( kCR_CRZVmT2_100  ,) \
  VALS( kCR_CRTopZjets   ,) \
  VALS( kCR_CRZXZjets    ,) \
  VALS( kCR_SSInc        ,) \
  VALS( kCR_PremT2       ,) \
  VALS( kCR_SRWHSS       ,) \
    /* DG These are the ones originally in SusyPlotter */ \
  VALS( kPR_NONE         ,) \
  VALS( kPR_SR6base      ,) \
  VALS( kPR_SR6          ,) \
  VALS( kPR_SR7base      ,) \
  VALS( kPR_SR7Nj        ,) \
  VALS( kPR_SR7NjZttVeto ,) \
  VALS( kPR_SR7NjPtTot   ,) \
  VALS( kPR_SR7NjMll     ,) \
  VALS( kPR_SR7          ,) \
  VALS( kPR_SR8base      ,) \
  VALS( kPR_CR8lpt       ,) \
    /* looser regions for fake control plots, same for SR9lpt */ \
  VALS( kPR_CR8ee        ,) \
  VALS( kPR_CR8mm        ,) \
  VALS( kPR_CR8mmMtww    ,) \
  VALS( kPR_CR8mmHt      ,) \
  VALS( kPR_SR8          ,) \
  VALS( kPR_SR9base      ,) \
  VALS( kPR_CR9lpt       ,) \
  VALS( kPR_SR9          ,) \
    /* for now only used in SusySelection::passSrSs */  \
  VALS( kWH_SRSS1        ,) \
  VALS( kWH_SRSS2        ,) \
  VALS( kWH_SRSS3        ,) \
  VALS( kWH_SRSS4        ,) \

  DECLARE_ENUM(SelectionRegions, SOME_ENUM)

} //end namespace susywh

#endif // end include guard
