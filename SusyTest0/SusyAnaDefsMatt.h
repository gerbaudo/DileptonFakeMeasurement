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
  CRSSInc,
  SR_WHSS,
  CR8lpt,
  CR8ee,
  CR8mm,
  CR8mmMtww,
  CR8mmHt,
};

#endif
