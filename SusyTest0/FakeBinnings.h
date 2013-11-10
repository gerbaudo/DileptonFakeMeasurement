// emacs -*- C++ -*-
#ifndef SUSY_FAKE_FAKEBINNINGS_H
#define SUSY_FAKE_FAKEBINNINGS_H



// binnings used to derive the fake estimate

namespace susy {
namespace fake {
const float FakePtbins[] = {10,20,35,100};
const int   nFakePtbins   = 3;
const float coarseFakePtbins[] = {10,20,35,100};
const int  nCoarseFakePtbins   = 3;
const float Etabins[] = {0,1.37,1.52,2.1,2.5};
const int  nEtabins = 4;
const float CoarseEtabins[] = {0,1.37,2.5};
const int  nCoarseEtabins = 2;
const float Metbins[] = {0, 20, 30, 40, 50, 60, 100, 200};
const int   nMetbins  = 7;
const float Jetbins[] = {-0.5, 0.5, 1.5, 2.5, 3.5};
const float nJetbins  = 4;
} // end namespace fake
} // end namespace susy

#endif
