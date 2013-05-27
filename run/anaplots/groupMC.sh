#!/bin/bash

# dropped all 'cfappend' samples (see Matt's original script for the complete list)
# /home/mrelich/workarea/SUSY2012/SusyAna2012/run/anaplots/groupMCwChargeFlip.sh

ending="_May27_n0139"
outdir="merged"
#addition="CFPtFix" 
addition=""
append=${ending}.AnaHists 
cfappend=${ending}.ChargeFlipHists


# data
hadd -f ${outdir}/data_physics_Egamma${append}${addition}.root \
    period?.physics_Egamma${append}.root
hadd -f ${outdir}/data_physics_Muons${append}${addition}.root \
    period?.physics_Muons${append}.root

# ttbar
hadd -f ${outdir}/top${append}${addition}.root \
    ttbar_LeptonFilter${append}.root \
    ttbarW${append}.root \
    ttbarZ${append}.root  \
    ttbarWj${append}.root \
    ttbarZj${append}.root  \
    SingleTopWtChanIncl${append}.root \

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

