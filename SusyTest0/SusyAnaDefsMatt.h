#ifndef SusyAnaDefsMatt_h
#define SusyAnaDefsMatt_h

////////////////////////////////////////////////////////
// Class to store common enums and binning for all of //
// the plots that I make for this analysis.           //
////////////////////////////////////////////////////////

//---------------------------------------------//
// Useful enums
//---------------------------------------------//

enum ControlRegion 
{
  // DG These are needed to compute the SF and the rates (pseudo t&p)
  CR_Real = 0,       // Real Z window
  CR_SideLow,        // Real Side band Lower
  CR_SideHigh,       // Real Side band higher
  CR_HF,             // HF b-bbar tag and probe
  CR_HF_high,        // HF b-bbar tag and probe with met < 80
  CR_LFZjet,         // LF Z+jet tag and probe
  CR_LFWjet,         // LF W+jet tag and probe
  CR_Conv,           // Zmumu tag and probe
  CR_CFlip,          // Charge Flip
  CR_MCHeavy,        // MC cr heavy flavor
  CR_MCLight,        // MC cr light flavor
  CR_MCConv,         // MC cr conversions
  CR_MCQCD,          // MC cr QCD = LF+HF
  CR_MCALL,          // MC All
  CR_MCReal,         // MC cr Real
  CR_MCNone,         // MC cr with no metrel cuts
  
  // DG These are the SR, which we need b/c we want to compute the compositions.
  CR_SRmT2a,
  CR_SRmT2b,
  CR_SRmT2c,
  CR_SRWWa,
  CR_SRWWb,
  CR_SRWWc,
  CR_SRZjets,
  
  CR_VRSS,
  CR_CRWWMet,
  CR_CRWWmT2,

  CR_CRTopMet,
  CR_CRTopmT2,
  CR_CRZVMet,
  CR_CRZVmT2_90, 
  CR_CRZVmT2_120, 
  CR_CRZVmT2_150, 
  CR_CRZVmT2_100, 

  CR_CRTopZjets,
  CR_CRZXZjets,

  // Added
  CR_SSInc,

  CR_PremT2,

  CR_SRWHSS,

  CR_N
};

static string CRNames[] = 
{
  "realCR",
  "realSideLow",
  "realSideHigh",
  "fakeHF",
  "fakeHF_high",
  "fakeLFZjet",
  "fakeLFWjet",
  "fakeConv",
  "fakeCFlip",
  "heavyMC",
  "lightMC",
  "convMC",
  "qcdMC",
  "allMC",
  "realMC",
  "noneMC",

  "SRmT2a",
  "SRmT2b",
  "SRmT2c",
  "SRWWa",
  "SRWWb",
  "SRWWc",
  "SRZjets",
  
  "VRSS",
  "CRWWMet",
  "CRWWmT2",

  "CRTopMet",
  "CRTopmT2",
  "CRZVMet",
  "CRZVmT2_90",
  "CRZVmT2_120",
  "CRZVmT2_150",
  "CRZVmT2_100",

  "CRTopZjets",
  "CRZXZjets",
  
  "CR_SSInc",
  
  "CRPremT2",
  "CR_WHSS"

};

static string CRLabels[] =
{
  "Z Tag and Probe",
  "Z Side Band: Low",
  "Z Side Band: High",
  "HF Tag and Probe",
  "HF Tag and Probe (high)",
  "LF Z+jet Tag and Probe",
  "LF W+jet Tag and Probe",
  "Z->#mu#mu Conversion",
  "Charge-Flip",
  "MC Truth Heavy",
  "MC Truth Light",
  "MC Truth Conversion",
  "MC Truth Jet Fakes",
  "MC Truth ALL",
  "MC Truth Real",
  "MC Truth",

  "Signal Region mT2a",
  "Signal Region mT2b",
  "Signal Region mT2c",
  "Signal Region WWa",
  "Signal Region WWb",
  "Signal Region WWc",

  "Validation Region SS",
  "Control Region WW Met",
  "Control Region WW mT2",
  "Control Region Top Met",
  "Control Region Top mT2",
  "Control Region ZV Met",
  "Control Region ZV mT2 90",
  "Control Region ZV mT2 120",
  "Control Region ZV mT2 150",
  "Control Region ZV mT2 100",

  "Control Region Top Z+jets",
  "Control Region Top ZX Z+jets",

  "SS Inclusive",

  "Pre-mT2 region",
  "Control Region WH SS"

};

    
// Lepton Types
enum LeptonType
{
  LT_EL = 0,        // Electron
  LT_MU,            // Muon
  LT_N              
};

static string LTNames[] = 
{
  "elec",
  "muon"
};

static string LTLabels[] =
{
  "Electron",
  "Muon"
};

// Dilepton Types for iterative corrections
enum DiLepPair
{
  DL_TT = 0,
  DL_TL,
  DL_LT,
  DL_LL,
  DL_ALL,
  DL_N
};

static string DLPNames[] =
{
  "TT",
  "TL",
  "LT",
  "LL",
  "ALL"
};

// Histogram types
enum HistType
{
  HT_DEN = 0,       // Denominator -- all
  HT_NUM,           // Numerator -- tight
  //HT_EFF,           // Efficiency -- Num/Den
  HT_N
};

static string HTNames[] = 
{
  "den",
  "num"
  //"eff"
};

// Separation into dilepton sign
enum DilepSign
{
  DS_SS = 0,
  DS_OS,
  DS_Both,
  DS_N
};

static string DSNames[] = 
{
  "SS",
  "OS",
  "Both"
};

// Signal regions

// DG These are the ones where you actually want to compute the fake
// prediction (i.e. your signal regions of the analysis, including the
// ones where you want to plot stuff to make checks.)

enum SignalRegion
{
  SRmT2a = 0,
  SRmT2b,
  SRmT2c,
  SRWWa,
  SRWWb,
  SRWWc,
  SRZjets,

  VRSS,
  CRWWMet,
  CRWWmT2,

  CRTopMet,
  CRTopmT2,
  CRZVMet,
  CRZVmT2_90, 
  CRZVmT2_120, 
  CRZVmT2_150, 
  CRZVmT2_100, 

  CRTopZjets,
  CRZXZjets,

  CRSSInc,
  
  CRPremT2,
  SR_WHSS,
  SR_N
};

static string SRNames[] =
{
  "SRmT2a",
  "SRmT2b",
  "SRmT2c",
  "SRWWa",
  "SRWWb",
  "SRWWc",
  "SRZjets",
  
  "VRSS",
  "CRWWMet",
  "CRWWmT2",

  "CRTopMet",
  "CRTopmT2",
  "CRZVMet",
  "CRZVmT2_90",
  "CRZVmT2_120",
  "CRZVmT2_150",
  "CRZVmT2_100",

  "CRTopZjets",
  "CRZXZjets",
  
  "CR_SSInc",

  "CRPremT2",
  "SR_WHSS"
  

};

static string SRProperNames[] =
{
  "Signal Region mT2a",
  "Signal Region mT2b",
  "Signal Region mT2c",
  "Signal Region WWa",
  "Signal Region WWb",
  "Signal Region WWc",
  "Signal Region Z+jets",

  "Validation Region SS",
  "Control Region WW Met",
  "Control Region WW mT2",
  "Control Region Top Met",
  "Control Region Top mT2",
  "Control Region ZV Met",
  "Control Region ZV mT2 90",
  "Control Region ZV mT2 120",
  "Control Region ZV mT2 150",
  "Control Region ZV mT2 100",

  "Control Region Top Z+jets",
  "Control Region Top ZX Z+jets",

  "Control Region SS Inclusive",
  
  "Pre-mT2 Region",
  "Signal region WH SS"
};

// Plotting regions
enum PlotRegion{

  PR_SRmT2a = 0,
  PR_SRmT2b,
  PR_SRmT2c,
  PR_SRWWa,
  PR_SRWWb,
  PR_SRWWc,
  PR_SRZjets,

  PR_VRSS,
  PR_CRWWMet,
  PR_CRWWmT2,

  PR_CRTopMet,
  PR_CRTopmT2a,
  PR_CRTopmT2b,
  PR_CRTopmT2c,
  PR_CRTopZjets,
  PR_CRZVMet,
  PR_CRZVmT2,

  PR_PreSRZjets,
  PR_CRZjets,

  PR_SRDavide,

  PR_N
};

static string PRNames[] = {
  "SRmT2a",
  "SRmT2b",
  "SRmT2c",
  "SRWWa",
  "SRWWb",
  "SRWWc",
  "SRZjets",

  "VRSS",
  "CRWWMet",
  "CRWWmT2",
  "CRTopMet",
  "CRTopmT2a",
  "CRTopmT2b",
  "CRTopmT2c",
  "CRTopZjets",
  "CRZVMet",
  "CRZVmT2",

  "PreSRZjets",
  "CRZjets",

  "SRDavide"
};

// Contamination or actual
enum RateInformation 
{
  RI_Actual,              // The actual rate
  RI_Contam,              // The contamination
  RI_All,                 // put both in there
  RI_N
};

static string RINames[] =
{
  "actual",
  "contam",
  "all"
};

// MC control regions
enum MCControlRegion
{
  MC_CR_RA = 0,         // Met rel < 30
  MC_CR_RB,             // 40 < Met rel < 100
  MC_CR_RC,             // OS Jet Veto Met rel < 30
  MC_CR_RD,             // OS Jet Veto 40 < Met rel < 100
  MC_CR_NONE,           // No met rel cut
  MC_CR_N
};

static string MCCRNames[] =
{
  "mcCRA",
  "mcCRB",
  "mcCRC",
  "mcCRD",
  "mcCRNONE"
};

static string MCCRLabels[] =
{
  "#slash{E}^{rel}_{T} < 30",
  "40 < #slash{E}^{rel}_{T} < 100",
  "OS Jet Veto #slash{E}^{rel}_{T} < 30",
  "OS Jet Veto 40 < #slash{E}^{rel}_{T} < 100",
  "No #slash{E}^{rel}_{T} Req"
};

// Lepton sources
enum LeptonSource
{
  LS_HF = 0,
  LS_LF,
  LS_Conv,
  LS_Real,
  LS_QCD,
  LS_Unk,
  LS_N
};

static string LSNames[] = 
{
  "heavy",
  "light",
  "conv",
  "real",
  "qcd",
  "unknown",
};

// lepton channel
enum Chan {
  Ch_all = 0,
  Ch_ee,
  Ch_mm,
  Ch_em,
  Ch_N
};

// Lepton chan names
static string chanNames[] = {
  "all",
  "ee",
  "mm",
  "em"
};


//---------------------------------------------//
// Binning for histograms
//---------------------------------------------//

const float Ptbins[]  = {0,10,15,20,25,30,50,70,100};
const int  nPtbins  = 8;

//const float FakePtbins[] = {0,10,15,20,25,30,35,100};
//const int  nFakePtbins   = 7;
//const float FakeSFPtbins[] = {10,15,20,30,40,100};
//const int  nFakeSFPtbins   = 5;
//const float RealPtbins[] = {10,15,20,25,30,35,100};
//const int  nRealPtbins   = 6;
const float RealPtbins[] = {10,15,20,25,30,35,45,55,70,100};
const int  nRealPtbins   = 9;
const float FakePtbins[] = {10,20,35,100};
const int  nFakePtbins   = 3;
//const float FakePtbins[] = {10,20,35,100};
//const int  nFakePtbins   = 3;
//const float FakePtbins[] = {10,15,20,35,50,100};
//const int  nFakePtbins   = 5;
//const float coarseFakePtbins[] = {10,20,35,100};
//const int  nCoarseFakePtbins   = 3;
const float lessCoarseFakePtbins[] = {10,15,20,35,100};
const int  nLessCoarseFakePtbins   = 4;
const float coarseFakePtbins[] = {10,20,35,100};
const int  nCoarseFakePtbins   = 3;
const float ConvFakePtbins[] = {10,15,20,25,100};
const int  nConvFakePtbins   = 4;
const float LFFakePtbins[] = {10,15,25,50,100};
const int   nLFFakePtbins = 4;
const float FineFakePtbins[] = {10,12,15,20,30,100};
const float nFineFakePtbins = 5;

const float MLLbins [] = {61,71,81,101,111,121};
const int  nMLLbins    = 5;

const float Etabins[] = {0,1.37,1.52,2.1,2.5};
const int  nEtabins = 4;

const float CoarseEtabins[] = {0,1.37,2.5};
const int  nCoarseEtabins = 2;

const float Typemin = -0.5;
const float Typemax = 22.5;
const int     nType = 23;

const float Originmin = -0.5;
const float Originmax = 42.5;
const int     nOrigin = 43;

const float Metbins[] = {0, 20, 30, 40, 50, 60, 100, 200};
const int   nMetbins  = 7;
const float coarseMetbins[] = {40, 100, 200};
const int   nCoarseMetbins  = 2;

const int nNPVbins = 30;
const float NPVmin = -0.5;
const float NPVmax = nNPVbins - 0.5;

const int nMUbins = 30;
const float MUmin = -0.5;
const float MUmax = nMUbins - 0.5;

const float Jetbins[] = {-0.5, 0.5, 1.5, 2.5, 3.5};
const float nJetbins  = 4;

const float OptIsomin = 0;
const float OptIsomax = 0.2;
const int nOptIsobins = 20;

const float OptSigmin = 0;
const float OptSigmax = 10;
const int nOptSigbins = 50;

const float OptZ0min = 0;
const float OptZ0max = 2;
const int nOptZ0bins = 40;


#endif
