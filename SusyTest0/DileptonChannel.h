// emacs -*- C++ -*-
#ifndef SUSY_WH_DILEPTONCHANNEL_H
#define SUSY_WH_DILEPTONCHANNEL_H

#include <string>

namespace susy {
namespace wh {

enum Chan                   { Ch_all = 0, Ch_ee, Ch_mm, Ch_em, Ch_N };
static std::string chanNames[] = {   "all",      "ee",  "mm",  "em"      };

}
}

#endif
