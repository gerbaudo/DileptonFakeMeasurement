// emacs -*- C++ -*-
#ifndef SUSY_FAKE_FAKELEPTONSOURCES_H
#define SUSY_FAKE_FAKELEPTONSOURCES_H

// Types and names of the possible sources of fake leptons

namespace susy {
namespace fake {

enum LeptonSource         { LS_HF = 0, LS_LF,   LS_Conv, LS_Real, LS_QCD, LS_Unk,  LS_N };
static string LSNames[] = { "heavy",   "light", "conv",  "real",  "qcd", "unknown"      };

}
}
#endif
