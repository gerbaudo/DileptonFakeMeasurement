// emacs -*- C++ -*-
#ifndef SUSY_WH_PLOTREGION_H
#define SUSY_WH_PLOTREGION_H

#include <string>
namespace susy {
namespace wh {

enum Region{
    PR_SR8base,
    PR_CR8lpt, PR_CR8ee, // looser regions for fake control plots, same for SR9lpt
    PR_CR8mm, PR_CR8mmMtww, PR_CR8mmHt,
    PR_SR8, PR_SR9base, PR_CR9lpt, PR_SR9,
    PR_SsEwk,
    PR_SsEwkLoose,
    CrZVfake1jee,
    CrZVfake2jee,
    CrZVfake1jem,
    CrZVfake2jem,
    Crfake1jem,
    Crfake2jem,
    CrZV1jmm,
    CrZV2jmm,
    Crfake1jmm,
    Crfake2jmm
};

const Region PlotRegions[] = {
    PR_SR8base,
    PR_CR8lpt, PR_CR8ee,
    PR_CR8mm, PR_CR8mmMtww, PR_CR8mmHt,
    PR_SR8, PR_SR9base, PR_CR9lpt, PR_SR9,
    PR_SsEwk, PR_SsEwkLoose,
    CrZVfake1jee,
    CrZVfake2jee,
    CrZVfake1jem,
    CrZVfake2jem,
    Crfake1jem,
    Crfake2jem,
    CrZV1jmm,
    CrZV2jmm,
    Crfake1jmm,
    Crfake2jmm
};
const size_t kNumberOfPlotRegions = sizeof(PlotRegions) / sizeof(PlotRegions[0]);
const string RegionNames[] =
{
  "sr8base"
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
};

inline std::string region2str(const Region &r) {return RegionNames[r];}

} // end namespace wh
} // end namespace susy
#endif
