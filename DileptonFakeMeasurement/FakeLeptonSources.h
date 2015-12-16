// emacs -*- C++ -*-
#ifndef SUSY_FAKE_FAKELEPTONSOURCES_H
#define SUSY_FAKE_FAKELEPTONSOURCES_H

#include <string>

namespace susy {
namespace fake {

enum LeptonSource              { LS_HF = 0, LS_LF,   LS_Conv, LS_Real, LS_QCD, LS_Unk,  LS_N };
static std::string LSNames[] = { "heavy",   "light", "conv",  "real",  "qcd", "unknown"      };

}
}
#endif
