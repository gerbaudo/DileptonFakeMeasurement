#include "SusyTest0/Systematics.h"

namespace susy {
namespace wh {


SusyNtSys sys2ntsys(const Systematic &s)
{
    SusyNtSys r = NtSys_N;
    switch(s) {
    case WH_CENTRAL          :  r =  NtSys_NOM         ; break;
    case WH_JESUP            :  r =  NtSys_JES_UP      ; break;
    case WH_JESDOWN          :  r =  NtSys_JES_DN      ; break;
    case WH_JER              :  r =  NtSys_JER         ; break;
    case WH_EESLOWUP         :  r =  NtSys_EES_LOW_UP  ; break;
    case WH_EESLOWDOWN       :  r =  NtSys_EES_LOW_DN  ; break;
    case WH_EESMATUP         :  r =  NtSys_EES_MAT_UP  ; break;
    case WH_EESMATDOWN       :  r =  NtSys_EES_MAT_DN  ; break;
    case WH_EESPSUP          :  r =  NtSys_EES_PS_UP   ; break;
    case WH_EESPSDOWN        :  r =  NtSys_EES_PS_DN   ; break;
    case WH_EERUP            :  r =  NtSys_EER_UP      ; break;
    case WH_EERDOWN          :  r =  NtSys_EER_DN      ; break;
    case WH_MESUP            :  r =  NtSys_MS_UP       ; break;
    case WH_MESDOWN          :  r =  NtSys_MS_DN       ; break;
    case WH_MIDUP            :  r =  NtSys_ID_UP       ; break;
    case WH_MIDDOWN          :  r =  NtSys_ID_DN       ; break;
    case WH_MEFFDOWN         :  r =  NtSys_EES_Z_DN    ; break;
    case WH_SCALESTUP        :  r =  NtSys_SCALEST_UP  ; break;
    case WH_SCALESTDOWN      :  r =  NtSys_SCALEST_DN  ; break;
    case WH_RESOST           :  r =  NtSys_RESOST      ; break;
    case WH_TESUP            :  r =  NtSys_TES_UP      ; break;
    case WH_TESDOWN          :  r =  NtSys_TES_DN      ; break;
    case WH_TTRIGSFUP        :  r =  NtSys_TRIGSF_EL_UP; break; // ? NtSys_TRIGSF_MU_UP
    case WH_TTRIGSFDOWN      :  r =  NtSys_TRIGSF_EL_DN; break; // ? NtSys_TRIGSF_MU_DN
    default : std::cout<<"cannot convert Systematic '"<<syst2str(s)<<"' to SusyNtSys; returning invalid NtSys_N"<<std::endl;
    }
    return r;
}
//-----------------------------------------
Systematic ntsys2sys(const SusyNtSys &s)
{
    int r=0;
    // here we assume that NtSys_N <= N(Systematic), which in general
    // is true because we have other syst (eg. fakes) in addition to
    // the SusyNt ones.
    while(r<NtSys_N) { if(s==sys2ntsys(static_cast<Systematic>(r))) break; else r++;}
    return static_cast<Systematic>(r);
}

} // wh
} // susy
