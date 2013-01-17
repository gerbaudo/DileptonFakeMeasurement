#!/bin/bash

# dropped all 'cfappend' samples (see Matt's original script for the complete list)
# /home/mrelich/workarea/SUSY2012/SusyAna2012/run/anaplots/groupMCwChargeFlip.sh

ending="_AE_Dec17_n0115"
ending="_Jan15_n0115"
outdir="merged"
#addition="CFPtFix" 
addition=""
append=${ending}.AnaHists 
cfappend=${ending}.ChargeFlipHists

# ttbar
#    PowhegPythia_AUET2BCT10_ttbar_LeptonFilter_AF2${append}.root \
hadd -f ${outdir}/top${append}${addition}.root \
    TtbarHadTauhad${append}.root \
    TtbarLeptHad${append}.root \
    TtbarLeptLept${append}.root \
    TtbarLeptTauhad${append}.root \
    TtbarLeptTaulept${append}.root \
    TtbarTauhadTauhad${append}.root \
    TtbarTauleptHad${append}.root \
    TtbarTauleptTauhad${append}.root \
    TtbarTauleptTaulept${append}.root \
    ttbarW${append}.root \
    ttbarZ${append}.root  \
    ttbarWj${append}.root \
    ttbarZj${append}.root  \
    SingleTopWtChanIncl${append}.root \


# Z+X
#hadd -f ${outdir}/ZX${append}${addition}.root \
#    Sherpa_CT10_Zee${append}.root \
#    Sherpa_CT10_Zmumu${append}.root \
#    ZeeNp0Excl_Mll10to60${append}.root \
#    ZeeNp1Excl_Mll10to60${append}.root \
#    ZeeNp2Excl_Mll10to60${append}.root \
#    ZeeNp3Excl_Mll10to60${append}.root \
#    ZeeNp4Excl_Mll10to60${append}.root \
#    ZmumuNp0Excl_Mll10to60${append}.root \
#    ZmumuNp1Excl_Mll10to60${append}.root \
#    ZmumuNp2Excl_Mll10to60${append}.root \
#    ZmumuNp3Excl_Mll10to60${append}.root \
#    ZmumuNp4Excl_Mll10to60${append}.root \
#    llll_ZZ${append}.root \
#    lllnu_WZ${append}.root \
#    llnunu_ZZ${append}.root \
#    ZZ4lep${append}.root \
#    ZZ4e${append}.root \
#    ZZ4mu${append}.root \
#    ZZ2e2mu${append}.root \
#    Sherpa_CT10_Zee${cfappend}.root \
#    Sherpa_CT10_Zmumu${cfappend}.root \
#    ZeeNp0Excl_Mll10to60${cfappend}.root \
#    ZeeNp1Excl_Mll10to60${cfappend}.root \
#    ZeeNp2Excl_Mll10to60${cfappend}.root \
#    ZeeNp3Excl_Mll10to60${cfappend}.root \
#    ZeeNp4Excl_Mll10to60${cfappend}.root \
#    ZmumuNp0Excl_Mll10to60${cfappend}.root \
#    ZmumuNp1Excl_Mll10to60${cfappend}.root \
#    ZmumuNp2Excl_Mll10to60${cfappend}.root \
#    ZmumuNp3Excl_Mll10to60${cfappend}.root \
#    ZmumuNp4Excl_Mll10to60${cfappend}.root \
#    llll_ZZ${cfappend}.root \
#    lllnu_WZ${cfappend}.root \
#    llnunu_ZZ${cfappend}.root \
#    ZZ4lep${cfappend}.root \
#    ZZ4e${cfappend}.root \
#    ZZ4mu${cfappend}.root \
#    ZZ2e2mu${cfappend}.root 


# Z+jet
hadd -f ${outdir}/Zjet${append}${addition}.root \
    Sherpa_CT10_Zee${append}.root \
    Sherpa_CT10_Zmumu${append}.root \
    ZeeNp0Excl_Mll10to60${append}.root \
    ZeeNp1Excl_Mll10to60${append}.root \
    ZeeNp2Excl_Mll10to60${append}.root \
    ZeeNp3Excl_Mll10to60${append}.root \
    ZeeNp4Excl_Mll10to60${append}.root \
    ZmumuNp0Excl_Mll10to60${append}.root \
    ZmumuNp1Excl_Mll10to60${append}.root \
    ZmumuNp2Excl_Mll10to60${append}.root \
    ZmumuNp3Excl_Mll10to60${append}.root \
    ZmumuNp4Excl_Mll10to60${append}.root \
    Sherpa_CT10_Ztautau${append}.root \
    ZtautauNp0Excl_Mll10to60${append}.root \
    ZtautauNp1Excl_Mll10to60${append}.root \
    ZtautauNp2Excl_Mll10to60${append}.root \
    ZtautauNp3Excl_Mll10to60${append}.root \


# Ztautau
#hadd -f ${outdir}/Ztautau${append}${addition}.root \
#    Sherpa_CT10_Ztautau${append}.root \
#    ZtautauNp0Excl_Mll10to60${append}.root \
#    ZtautauNp1Excl_Mll10to60${append}.root \
#    ZtautauNp2Excl_Mll10to60${append}.root \
#    ZtautauNp3Excl_Mll10to60${append}.root \
#    Sherpa_CT10_Ztautau${cfappend}.root \
#    ZtautauNp0Excl_Mll10to60${cfappend}.root \
#    ZtautauNp1Excl_Mll10to60${cfappend}.root \
#    ZtautauNp2Excl_Mll10to60${cfappend}.root \
#    ZtautauNp3Excl_Mll10to60${cfappend}.root 

 
# Diboson (Z)
hadd -f ${outdir}/ZZ${append}${addition}.root \
    llll_ZZ${append}.root \
    llnunu_ZZ${append}.root \
    ZZ4lep${append}.root \
    ZZ4e${append}.root \
    ZZ4mu${append}.root \
    ZZ2e2mu${append}.root \


hadd -f ${outdir}/WZ${append}${addition}.root \
    lllnu_WZ${append}.root \


# WW
hadd -f ${outdir}/WW${append}${addition}.root \
    llnunu_WW${append}.root \
    llnunu_SS_EW6${append}.root \
    llnunujj_SS${append}.root \
    WpWmenuenu${append}.root \
    WpWmenumunu${append}.root \
    WpWmenutaunu${append}.root \
    WpWmmunumunu${append}.root \
    WpWmmunutaunu${append}.root \
    WpWmtaunuenu${append}.root \
    WpWmtaunumunu${append}.root \
    WpWmtaunutaunu${append}.root \


#hadd -f ${outdir}/HF${append}${addition}.root \
#    bbTomu20${append}.root \
#    ccTomu20${append}.root \
#    bbTomu20${cfappend}.root \
#    ccTomu20${cfappend}.root 

#WpWmmunuenu${append}.root <-- Doesn't exist yet

    #enugammaPt10${append}.root 
    #munugammaPt10${append}.root
    #taunugammaPt10${append}.root
    #eegammaPt10${append}.root 
    #mumugammaPt10${append}.root 
    #tautaugammaPt10${append}.root 

#hadd -f ${outdir}/totalMC${append}${addition}.root \
#    top${append}${addition}.root \
#    WW${append}${addition}.root \
#    Zdib${append}${addition}.root \
#    Ztautau${append}${addition}.root \
#    Zjet${append}${addition}.root \
#    data_Dec4_n0111.FakeHists.root
