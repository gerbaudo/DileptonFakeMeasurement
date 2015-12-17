// emacs -*- C++ -*-
#ifndef SUSY_FAKE_FAKEREGIONS_H
#define SUSY_FAKE_FAKEREGIONS_H

#include <string>


// binnings used to derive the fake estimate

namespace susy {
namespace fake {

enum Region
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
  CR_SSInc1j,
  CR_SRWHSS,
  CR_CR8lpt,
  CR_CR8ee,
  CR_CR8mm,
  CR_CR8mmMtww,
  CR_CR8mmHt,
  CR_CR9lpt,
  CR_SsEwk,
  CR_SsEwkLoose,
  CR_SsEwkLea,
  CR_WHZVfake1jee,
  CR_WHZVfake2jee,
  CR_WHZVfake1jem,
  CR_WHZVfake2jem,
  CR_WHfake1jem,
  CR_WHfake2jem,
  CR_WHZV1jmm,
  CR_WHZV2jmm,
  CR_WHfake1jmm,
  CR_WHfake2jmm,

  CR_WHZVfake1j,
  CR_WHZVfake2j,
  CR_WHfake1j,
  CR_WHfake2j,
  CR_WHZV1j,
  CR_WHZV2j,

  CR_SRWH1j,
  CR_SRWH2j,
  CR_SRWHnoMlj, // same as above but without n-jet and mljj requirements
  CR_HF_SS, // see MeasureFakeRate2::passHFCR_Ss
  CR_HF_mme, // see MeasureFakeRate2::passZvetoPlusJetsCR
  CR_emu, // for hlfv  
  CR_razor0j, // for Danny's study
  CR_razor1j,
};

const Region SignalRegions[] = {
    CR_SSInc1j,
    CR_SRWHnoMlj
/*
    CR_SSInc, CR_SSInc1j, CR_SRWHSS,
    CR_CR8lpt, CR_CR8ee, CR_CR8mm, CR_CR8mmMtww, CR_CR8mmHt,
    CR_CR9lpt,
    CR_SsEwk, CR_SsEwkLoose, CR_SsEwkLea,
    CR_WHZVfake1jee,
    CR_WHZVfake2jee,
    CR_WHZVfake1jem,
    CR_WHZVfake2jem,
    CR_WHfake1jem,
    CR_WHfake2jem,
    CR_WHZV1jmm,
    CR_WHZV2jmm,
    CR_WHfake1jmm,
    CR_WHfake2jmm,

    CR_WHZVfake1j,
    CR_WHZVfake2j,
    CR_WHfake1j,
    CR_WHfake2j,
    CR_WHZV1j,
    CR_WHZV2j,

    CR_SRWH1j,
    CR_SRWH2j,
    CR_SRWHnoMlj
*/
};
const int NumberOfSignalRegions = sizeof(SignalRegions) / sizeof(SignalRegions[0]);
const std::string RegionNames[] =
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
  "CR_SSInc1j",
  "CR_WHSS",
  "CR_CR8lpt",
  "CR_CR8ee",
  "CR_CR8mm",
  "CR_CR8mmMtww",
  "CR_CR8mmHt",
  "CR_CR9lpt",
  "CR_SsEwk",
  "CR_SsEwkLoose",
  "CR_SsEwkLea",
  "CR_WHZVfake1jee",
  "CR_WHZVfake2jee",
  "CR_WHZVfake1jem",
  "CR_WHZVfake2jem",
  "CR_WHfake1jem",
  "CR_WHfake2jem",
  "CR_WHZV1jmm",
  "CR_WHZV2jmm",
  "CR_WHfake1jmm",
  "CR_WHfake2jmm",

  "CR_WHZVfake1j",
  "CR_WHZVfake2j",
  "CR_WHfake1j",
  "CR_WHfake2j",
  "CR_WHZV1j",
  "CR_WHZV2j",

  "CR_SRWH1j",
  "CR_SRWH2j",
  "CR_SRWHnoMlj",
  "CR_HF_SS",
  "CR_HF_mme",
  "CR_emu",
  "CR_razor0j",
  "CR_razor1j"
};

inline std::string region2str(const Region &r) {return RegionNames[r];}

} // end namespace fake
} // end namespace susy

#endif
