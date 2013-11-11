#ifndef SusyAnaDefsMatt_h
#define SusyAnaDefsMatt_h

////////////////////////////////////////////////////////
// Class to store common enums and binning for all of //
// the plots that I make for this analysis.           //
////////////////////////////////////////////////////////

//---------------------------------------------//
// Useful enums
//---------------------------------------------//

// Lepton Types
enum LeptonType            { LT_EL = 0,  LT_MU, LT_N };
static string LTNames[] =  { "elec",     "muon"      };
// Dilepton Types for iterative corrections
enum DiLepPair             { DL_TT = 0, DL_TL, DL_LT, DL_LL, DL_ALL, DL_N };

// Signal regions; these were used in FinalNewFake to determine the
// signal regions where we want to compute the fake prediction.
// Now we are only parsing the enum in buildWeightedMatrix.py.
// There is no need for an additional enum, and we will drop this.
enum SignalRegion
{
  CR_SSInc,
  CR_WHSS,
  CR_CR8lpt,
  CR_CR8ee,
  CR_CR8mm,
  CR_CR8mmMtww,
  CR_CR8mmHt,
  CR_SsEwk,
  CR_SsEwkLoose,
};

#endif
