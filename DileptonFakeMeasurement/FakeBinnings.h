// emacs -*- C++ -*-
#ifndef SUSY_FAKE_FAKEBINNINGS_H
#define SUSY_FAKE_FAKEBINNINGS_H

// binnings used to derive the fake estimate

namespace susy {
namespace fake {
const float FakePtbins[]       = {10,20,35,100};
const float coarseFakePtbins[] = {10,20,35,100};
const float Etabins[]          = {0,1.37,1.52,2.1,2.5};
const float CoarseEtabins[]    = {0,1.37,2.5};
const float Metbins[]          = {0, 20, 30, 40, 50, 60, 100, 200};
const float Jetbins[]          = {-0.5, 0.5, 1.5, 2.5, 3.5};
const float flavBins[]         = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
const int   nFakePtbins        = sizeof(FakePtbins)       / sizeof(FakePtbins[0])       - 1;
const int   nCoarseFakePtbins  = sizeof(coarseFakePtbins) / sizeof(coarseFakePtbins[0]) - 1;
const int   nEtabins           = sizeof(Etabins)          / sizeof(Etabins[0])          - 1;
const int   nCoarseEtabins     = sizeof(CoarseEtabins)    / sizeof(CoarseEtabins[0])    - 1;
const int   nMetbins           = sizeof(Metbins)          / sizeof(Metbins[0])          - 1;
const int   nJetbins           = sizeof(Jetbins)          / sizeof(Jetbins[0])          - 1;
const int   nFlavBins          = sizeof(flavBins)         / sizeof(flavBins[0])         - 1;
} // end namespace fake
} // end namespace susy

#endif
