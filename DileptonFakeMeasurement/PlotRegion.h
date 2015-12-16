// emacs -*- C++ -*-
#ifndef SUSY_WH_PLOTREGION_H
#define SUSY_WH_PLOTREGION_H

#include <string>
namespace susy {
namespace wh {

enum Region{
    PR_SR8base,
    PR_CRSsInc1j,
    PR_CR8lpt, PR_CR8ee, // looser regions for fake control plots, same for SR9lpt
    PR_CR8mm, PR_CR8mmMtww, PR_CR8mmHt,
    PR_SR8, PR_SR9base, PR_CR9lpt, PR_SR9,
    PR_SsEwk,
    PR_SsEwkLoose,
    PR_SsEwkLea,
    CrZVfake1jee,
    CrZVfake2jee,
    CrZVfake1jem,
    CrZVfake2jem,
    Crfake1jem,
    Crfake2jem,
    CrZV1jmm,
    CrZV2jmm,
    Crfake1jmm,
    Crfake2jmm,

    CrZVfake1j,
    CrZVfake2j,
    Crfake1j,
    Crfake2j,
    CrZV1j,
    CrZV2j,

    SrWh1j,
    SrWh2j

};

const Region PlotRegions[] = {
    PR_SR8base,
    PR_CRSsInc1j,
    PR_CR8lpt, PR_CR8ee,
    PR_CR8mm, PR_CR8mmMtww, PR_CR8mmHt,
    PR_SR8, PR_SR9base, PR_CR9lpt, PR_SR9,
    PR_SsEwk, PR_SsEwkLoose, PR_SsEwkLea,
    CrZVfake1jee,
    CrZVfake2jee,
    CrZVfake1jem,
    CrZVfake2jem,
    Crfake1jem,
    Crfake2jem,
    CrZV1jmm,
    CrZV2jmm,
    Crfake1jmm,
    Crfake2jmm,

    CrZVfake1j,
    CrZVfake2j,
    Crfake1j,
    Crfake2j,
    CrZV1j,
    CrZV2j,
    SrWh1j,
    SrWh2j
};
const size_t kNumberOfPlotRegions = sizeof(PlotRegions) / sizeof(PlotRegions[0]);
const string RegionNames[] =
{
  "sr8base"
  ,"crSsInc1j"
  ,"cr8lpt"
  ,"cr8lptee"
  ,"cr8lptmm"
  ,"cr8lptmmMtww"
  ,"cr8lptmmHt"
  ,"sr8"
  ,"sr9base"
  ,"cr9lpt"
  ,"sr9"
  ,"srSsEwk"
  ,"crSsEwkLoose"
  ,"crSsEwkLea"
  ,"crZVfake1jee"
  ,"crZVfake2jee"
  ,"crZVfake1jem"
  ,"crZVfake2jem"
  ,"crfake1jem"
  ,"crfake2jem"
  ,"crZV1jmm"
  ,"crZV2jmm"
  ,"crfake1jmm"
  ,"crfake2jmm"

  ,"crZVfake1j"
  ,"crZVfake2j"
  ,"crfake1j"
  ,"crfake2j"
  ,"crZV1j"
  ,"crZV2j"
  ,"srWh1j"
  ,"srWh2j"
};

inline std::string region2str(const Region &r) {return RegionNames[r];}

} // end namespace wh
} // end namespace susy
#endif
