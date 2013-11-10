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
// Lepton sources
enum LeptonSource         { LS_HF = 0, LS_LF,   LS_Conv, LS_Real, LS_QCD, LS_Unk,  LS_N };
static string LSNames[] = { "heavy",   "light", "conv",  "real",  "qcd", "unknown"      };
// lepton channel
enum Chan                   { Ch_all = 0, Ch_ee, Ch_mm, Ch_em, Ch_N };
static string chanNames[] = {   "all",      "ee",  "mm",  "em"      };

enum ControlRegion // these are the regions used in MeasureFakeRate2
{
  // DG These are needed to compute the SF and the rates (pseudo t&p)
  CR_Real = 0,       // Real Z window
  CR_SideLow,        // Real Side band Lower
  CR_SideHigh,       // Real Side band higher
  CR_HF,             // HF b-bbar tag and probe
  CR_HF_high,        // HF b-bbar tag and probe with met < 80
  CR_Conv,           // Zmumu tag and probe
  CR_MCConv,         // MC cr conversions
  CR_MCQCD,          // MC cr QCD = LF+HF
  CR_MCReal,         // MC cr Real
  // DG These are the SR, which we need b/c we want to compute the compositions.
  CR_SSInc,
  CR_SRWHSS,
  CR_CR8lpt,
  CR_CR8ee,
  CR_CR8mm,
  CR_CR8mmMtww,
  CR_CR8mmHt,
  //  CR_N
};

static string CRNames[] =
{
  "realCR",
  "realSideLow",
  "realSideHigh",
  "fakeHF",
  "fakeHF_high",
  "fakeConv",
  "convMC",
  "qcdMC",
  "realMC",
  "CR_SSInc",
  "CR_WHSS",
  "CR_CR8lpt",
  "CR_CR8ee",
  "CR_CR8mm",
  "CR_CR8mmMtww",
  "CR_CR8mmHt",

};
// Signal regions; these are used in FinalNewFake
// DG These are the ones where you actually want to compute the fake
// prediction (i.e. your signal regions of the analysis, including the
// ones where you want to plot stuff to make checks.)
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
