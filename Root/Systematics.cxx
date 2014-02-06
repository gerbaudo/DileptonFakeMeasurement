#include "SusyTest0/Systematics.h"

#include <cassert>

namespace susy {
namespace wh {

//-----------------------------------------
Systematic ntsys2sys(const SusyNtSys &s)
{
    Systematic r = WH_CENTRAL;
    switch(s) {
    case NtSys_NOM             :  r =  WH_CENTRAL      ; break;
    case NtSys_EES_Z_UP        :  r =  WH_EESZUP       ; break;
    case NtSys_EES_Z_DN        :  r =  WH_EESZDOWN     ; break;
    case NtSys_EES_MAT_UP      :  r =  WH_EESMATUP     ; break;
    case NtSys_EES_MAT_DN      :  r =  WH_EESMATDOWN   ; break;
    case NtSys_EES_PS_UP       :  r =  WH_EESPSUP      ; break;
    case NtSys_EES_PS_DN       :  r =  WH_EESPSDOWN    ; break; 
    case NtSys_EES_LOW_UP      :  r =  WH_EESLOWUP     ; break;
    case NtSys_EES_LOW_DN      :  r =  WH_EESLOWDOWN   ; break;
    case NtSys_EER_UP          :  r =  WH_EERUP        ; break;
    case NtSys_EER_DN          :  r =  WH_EERDOWN      ; break; 
    case NtSys_MS_UP           :  r =  WH_MESUP        ; break;
    case NtSys_MS_DN           :  r =  WH_MESDOWN      ; break;
    case NtSys_ID_UP           :  r =  WH_MIDUP        ; break;
    case NtSys_ID_DN           :  r =  WH_MIDDOWN      ; break;
    case NtSys_JES_UP          :  r =  WH_JESUP        ; break;
    case NtSys_JES_DN          :  r =  WH_JESDOWN      ; break;
    case NtSys_JER             :  r =  WH_JER          ; break;
    case NtSys_SCALEST_UP      :  r =  WH_SCALESTUP    ; break;
    case NtSys_SCALEST_DN      :  r =  WH_SCALESTDOWN  ; break;
    case NtSys_RESOST          :  r =  WH_RESOST       ; break;
    case NtSys_TRIGSF_EL_UP    :  r =  WH_TTRIGSFUP    ; break;
    case NtSys_TRIGSF_EL_DN    :  r =  WH_TTRIGSFDOWN  ; break;
    case NtSys_TRIGSF_MU_UP    :  /* undefined ?? */   ; break;
    case NtSys_TRIGSF_MU_DN    :  /* undefined ?? */   ; break;
    case NtSys_TES_UP          :  r =  WH_TESUP        ; break;
    case NtSys_TES_DN          :  r =  WH_TESDOWN      ; break;
    case NtSys_N               : assert(false)         ; break; // perhaps throw an exception instead
        // no default, so that the compiler will warn us of un-handled cases
    }
    return r;
}
//-----------------------------------------
SusyNtSys sys2ntsys(const Systematic &s)
{
    int r=0;
    // Here we assume that NtSys_N <= N(Systematic); in general this
    // is true because we have other syst (eg. fakes) in addition to
    // the SusyNt ones.
    while(r<NtSys_N) { if(s==ntsys2sys(static_cast<SusyNtSys>(r))) break; else r++;}
    return static_cast<SusyNtSys>(r);
}
//-----------------------------------------

} // wh
} // susy
