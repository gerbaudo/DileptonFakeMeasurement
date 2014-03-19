// emacs -*- C++ -*-
#ifndef SUSY_WH_SYSTEMATICS_H
#define SUSY_WH_SYSTEMATICS_H

/*

List of possible systematic uncertainties for the WH analysis
  
Essentially a copy of HistFitterTree/Systematics_naming.txt, with a
few additional helper functions to handle names and to convert between
susy::wh::Systematic <-> SusyNtSys

davide.gerbaudo@gmail.com
Feb 2014
*/

#include "SusyNtuple/SusyDefs.h"
#include <string>

namespace susy {
namespace wh {

enum Systematic {
    WH_CENTRAL             // Central value
    ,WH_STAT               // Statistical Uncertainty (no tree/weight needed for Limit setting)
    ,WH_JESUP              // Positive Jet energy scale
    ,WH_JESDOWN            // Negative Jet energy scale
    ,WH_JER                // Jet energy resolution
    ,WH_EESLOWUP           // Positive shift in electron energy scale (LOW)
    ,WH_EESLOWDOWN         // Negative shift in electron energy scale (LOW)
    ,WH_EESMATUP           // Positive shift in electron energy scale (Mat)
    ,WH_EESMATDOWN         // Negative shift in electron energy scale (Mat)
    ,WH_EESPSUP            // Positive shift in electron energy scale (PS)
    ,WH_EESPSDOWN          // Negative shift in electron energy scale (PS)
    ,WH_EESZUP             // Positive shift in electron energy scale (Zee)
    ,WH_EESZDOWN           // Negative shift in electron energy scale (Zee)
    ,WH_EERUP              // Positive shift in electron energy resolution
    ,WH_EERDOWN            // Negative shift in electron energy resolution
    ,WH_MESUP              // Positive shift in muon energy scale
    ,WH_MESDOWN            // Negative shift in muon energy scale
    ,WH_MIDUP              // Positive shift in muon ID resolution
    ,WH_MIDDOWN            // Negative shift in muon ID resolution
    ,WH_MMSUP              // Positive shift in muon MS resolution
    ,WH_MMSDOWN            // Negative shift in muon MS resolution
    ,WH_ESFUP              // Positive shift in electron efficiency
    ,WH_ESFDOWN            // Negative shift in electron efficiency
    ,WH_MEFFUP             // Positive shift in muon efficiency
    ,WH_MEFFDOWN           // Negative shift in muon efficiency
    ,WH_ETRIGREWUP         // Positive shift in electron trigger weights
    ,WH_ETRIGREWDOWN       // Negative shift in electron trigger weights
    ,WH_MTRIGREWUP         // Positive shift in muon trigger weights
    ,WH_MTRIGREWDOWN       // Negative shift in muon trigger weights
    ,WH_BJETUP             // Positive shift in btag scale factor
    ,WH_BJETDOWN           // Negative shift in btag scale factor
    ,WH_CJETUP             // Positive shift in ctag scale factor
    ,WH_CJETDOWN           // Negative shift in ctag scale factor
    ,WH_BMISTAGUP          // Positive shift in ltag scale factor
    ,WH_BMISTAGDOWN        // Negative shift in ltag scale factor
    ,WH_SCALESTUP          // Positive shift in MET soft term scale
    ,WH_SCALESTDOWN        // Negative shift in MET soft term scale
    ,WH_RESOST             // MET soft term resolution 
    ,WH_GEN                // Uncertainty due to generator (Alpgen versus Sherpa,...) and patronshower/hadronization (Herwig versus Pythia) (possibly includes PDF choice as well, i.e. CTEQ versus HERA,...) 
    ,WH_GENUP              // Positive shift due to generator parameter variations (ktfac, qfac, ptmin, Iqopt, scale); includes also shift in ISR/FSR
    ,WH_GENDOWN            // Negative shift due to generator parameter variations (ktfac, qfac, ptmin, Iqopt, scale) ;includes also shift in ISR/FSR
    ,WH_PDFERRUP           // Positive shift due to PDF errorset
    ,WH_PDFERRDOWN         // Negative shift due to PDF errorset
    ,WH_BKGMETHODUP        // Positive shift due to background estimation method (chargeflip, jetveto, fakeweights,...) 
    ,WH_BKGMETHODDOWN      // Negative shift due to background estimation method (chargeflip, jetveto, fakeweights,...) 
    ,WH_XSUP               // Positive shift in theoretical cross-section uncertainty
    ,WH_XSDOWN             // Negative shift in theoretical cross-section uncertainty
    ,WH_LUMI               // Uncertainty due to luminosity (no tree/weight needed for Limit setting)
    ,WH_PILEUPUP           // Positive shift for mu 
    ,WH_PILEUPDOWN         // Negative shift for mu
    ,WH_JVFUP              // Positive shift in Jet Vertex Fraction
    ,WH_JVFDOWN            // Negative shift in Jet Vertex Fraction
    ,WH_TESUP              // Positive shift in TES variation                                             
    ,WH_TESDOWN            // Negative shift in TES variation
    ,WH_TIDSFUP            // Positive shift in tauID scale factor                                      
    ,WH_TIDSFDOWN          // Negative shift in tauID scale factor
    ,WH_TEVSFUP            // Positive shift in tauEVeto scale factor                                   
    ,WH_TEVSFDOWN          // Negative shift in tauEVeto scale factor
    ,WH_TTRIGSFUP          // Positive shift in tau trigger scale factor                              
    ,WH_TTRIGSFDOWN        // Negative shift in tau trigger scale factor
    ,WH_TFAKESFUP          // Positive shift in fake tau scale factors                                
    ,WH_TFAKESFDOWN        // Negative shift in fake tau scale factors  
    ,WH_FakeTauBGSyst      // Uncertainty on fake tau background estimation
};

const std::string SystematicNames[] = {
    "WH_CENTRAL"
    ,"WH_STAT"
    ,"WH_JESUP"
    ,"WH_JESDOWN"
    ,"WH_JER"
    ,"WH_EESLOWUP"
    ,"WH_EESLOWDOWN"
    ,"WH_EESMATUP"
    ,"WH_EESMATDOWN"
    ,"WH_EESPSUP"
    ,"WH_EESPSDOWN"
    ,"WH_EESZUP"
    ,"WH_EESZDOWN"
    ,"WH_EERUP"
    ,"WH_EERDOWN"
    ,"WH_MESUP"
    ,"WH_MESDOWN"
    ,"WH_MIDUP"
    ,"WH_MIDDOWN"
    ,"WH_MMSUP"
    ,"WH_MMSDOWN"
    ,"WH_ESFUP"
    ,"WH_ESFDOWN"
    ,"WH_MEFFUP"
    ,"WH_MEFFDOWN"
    ,"WH_ETRIGREWUP"
    ,"WH_ETRIGREWDOWN"
    ,"WH_MTRIGREWUP"
    ,"WH_MTRIGREWDOWN"
    ,"WH_BJETUP"
    ,"WH_BJETDOWN"
    ,"WH_CJETUP"
    ,"WH_CJETDOWN"
    ,"WH_BMISTAGUP"
    ,"WH_BMISTAGDOWN"
    ,"WH_SCALESTUP"
    ,"WH_SCALESTDOWN"
    ,"WH_RESOST"
    ,"WH_GEN"
    ,"WH_GENUP"
    ,"WH_GENDOWN"
    ,"WH_PDFERRUP"
    ,"WH_PDFERRDOWN"
    ,"WH_BKGMETHODUP"
    ,"WH_BKGMETHODDOWN"
    ,"WH_XSUP"
    ,"WH_XSDOWN"
    ,"WH_LUMI"
    ,"WH_PILEUPUP"
    ,"WH_PILEUPDOWN"
    ,"WH_JVFUP"
    ,"WH_JVFDOWN"
    ,"WH_TESUP"
    ,"WH_TESDOWN"
    ,"WH_TIDSFUP"
    ,"WH_TIDSFDOWN"
    ,"WH_TEVSFUP"
    ,"WH_TEVSFDOWN"
    ,"WH_TTRIGSFUP"
    ,"WH_TTRIGSFDOWN"
    ,"WH_TFAKESFUP"
    ,"WH_TFAKESFDOWN"
    ,"WH_FakeTauBGSyst"
};


inline bool isValid(const Systematic &s) { return s>=WH_CENTRAL && s<=WH_FakeTauBGSyst; }
inline bool isValid(const SusyNtSys &s)  { return s>=NtSys_NOM  && s< NtSys_N; }
inline std::string syst2str(const Systematic &s) { return isValid(s) ? SystematicNames[s] : "unknown"; }
inline std::string syst2str(const SusyNtSys &s)  { return isValid(s) ? SusyNtSystNames[s] : "unknown"; }
SusyNtSys sys2ntsys(const Systematic &s);
Systematic ntsys2sys(const SusyNtSys &s);
BTagSys sys2ntbsys(const Systematic &s); //!< convert to SusyDef::BTagSys; return nominal if not a btag-related sys
} // wh
} // susy
#endif
