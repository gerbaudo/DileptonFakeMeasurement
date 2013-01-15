#ifndef SusyAnaDefs_h
#define SusyAnaDefs_h

////////////////////////////////////////////////////////
// Class to store common enums and binning for all of //
// the plots that I make for this analysis.           //
////////////////////////////////////////////////////////

//---------------------------------------------//
// Useful enums
//---------------------------------------------//

// Control Regions
enum ControlRegion
{
  CR_RA = 0,        // Real CRA: Z window
  CR_RB,            // Real CRA: just left of Z window
  CR_RC,            // Real CRA: just right of Z window
  CR_RD,            // Real CRB: Z window for Opposite Flavor
  CR_RE,            // Real CRC: Z window for Same-sign
  CR_FA,            // Fake CRA: Single Lepton
  CR_FB,            // Fake CRB: SS Dilepton
  CR_FC_none,       // Fake CRC: Tag in jet, probe away from jet, No met rel check
  CR_FC_low,        // Fake CRC: Tag in jet, probe away from jet, MetRel < 30
  CR_FC_int,        // Fake CRC: Tag in jet, probe away from jet, 40<MetRel<100
  CR_FD,            // Fake CRD: Photon+jet
  CR_FE,            // Fake CRE: Z->mumu+photon
  CR_FE_true,       // Fake CRE: Z->mumu+photon requiring true lepton in MC
  CR_FF,            // Fake CRF: Z+jets for LF
  CR_FF_true,       // Fake CRF: Z+jets for LF requireing true lepton in MC
  CR_FGignac,       // Fake CR to match Matt Gignac's HF
  CR_FG,            // Fake CR to match Matt Gignac's LF EXACTLY
  CR_FH,            // Attempt at W+jet CR
  CR_N
};
   
static string CRNames[] =
{
  "realCRA",
  "realCRB",
  "realCRC",
  "realCRD",
  "realCRE",
  "fakeCRA",
  "fakeCRB",
  "fakeCRC_none",
  "fakeCRC_low",
  "fakeCRC_int",
  "fakeCRD",
  "fakeCRE",
  "fakeCRE_true",
  "fakeCRF",
  "fakeCRF_true",
  "fakeCRGignac",
  "fakeCRG",
  "fakeCRH"
};

static string CRLabels[] =
{
  "Z tag and probe",
  "Lower Z sideband",
  "Upper Z sideband",
  "Z tag and probe for OF",
  "Single Lepton CR",
  "SS Dilepton CR",
  "HF Tag and Probe",
  "HF Tag and Probe #slash{E}^{rel}_{T} < 30",
  "HF Tag and Probe 40 < #slash{E}^{rel}_{T} < 100",  
  "#gamma+jet CR",
  "Conversion CR",
  "Conversion CR True"
  "Light Flavor CR",
  "Light Flavor CR True",
  "Matt Gignac's HF",
  "Matt Gignac's LF",
  "W+jet LF"
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
enum SignalRegion
{
  SR1 = 0,   
  SR2, 
  SR3,
  SR4,
  SR4b,
  SR5,
  SRNONE,
  VR1,  
  VR2,
  VR3,
  WWCR1,
  WWCR2,
  WWCR3,
  TOPCR,
  BR1,
  BR2,
  BR3,
  BR4,
  SR_N
};

static string SRNames[] =
{
  "SR1",
  "SR2",
  "SR3",
  "SR4",
  "SR4b",
  "SR5",
  "SRNONE",
  "VR1",
  "VR2",
  "VR3",
  "WWCR1",
  "WWCR2",
  "WWCR3",
  "TOPCR",
  "BR1",
  "BR2",
  "BR3",
  "BR4"
};

static string SRProperNames[] =
{
  "OS Jet Veto",
  "SS Jet Veto",
  "2 Jet",
  "m_{T2}",
  "m_{T2} b",
  "SR5",
  "SRNONE",
  "Validation Region 1",
  "Validation Region 2",
  "Validation Region 3",
  "WW Control Region 1",
  "WW Control Region 2",
  "WW Control Region 3",
  "Top Control Region",
  "Bonus Region 1",
  "Bonus Region 2"
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
  LS_QCD,
  LS_Conv,
  LS_Real,
  LS_N
};

static string LSNames[] = 
{
  "heavy",
  "light",
  "qcd",
  "conv",
  "real"
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
const float FakePtbins[] = {10,15,20,35,100};
const int  nFakePtbins   = 4;
const float coarseFakePtbins[] = {10,20,30,100};
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

const float Metbins[] = {40, 50, 60, 100, 200};
const int   nMetbins  = 4;

const int nNPVbins = 30;
const float NPVmin = -0.5;
const float NPVmax = nNPVbins - 0.5;

#endif
