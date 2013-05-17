#!/bin/bash


append="_AE_Dec17_n0115"
append="_May17_n0139"


# only one option supported: --1fb to use the 1/fb weights
#option="--1fb"
#option="--AD"

echo "Appending " ${append} " Run option: " ${option}

DATA=(

    # Data
    periodA_Egamma periodA_Muons
    periodB_Egamma periodB_Muons
    periodC_Egamma periodC_Muons
    periodD_Egamma periodD_Muons
    periodE_Egamma periodE_Muons
    periodG_Egamma periodG_Muons
    periodH_Egamma periodH_Muons
    periodI_Egamma periodI_Muons
    periodJ_Egamma periodJ_Muons
    periodL_Egamma periodL_Muons
)

MC=(
    # New Z+Jets
    #ZeeLightJets_AF2 ZmumuLightJets_AF2 ZtautauLightJets_AF2
    #ZeeHeavyJets_AF2 ZmumuHeavyJets_AF2 ZtautauHeavyJets_AF2

    # Core Samples
    Sherpa_CT10_Zee
    Sherpa_CT10_Zmumu
    Sherpa_CT10_Ztautau

    # Z+jets
    #ZeeNp0 ZeeNp1 ZeeNp2 ZeeNp3 ZeeNp4 Zweep5
    #ZmumuNp0 ZmumuNp1 ZmumuNp2 ZmumuNp3 ZmumuNp4 ZmumuNp5
    #ZtautauNp0 ZtautauNp1 ZtautauNp2 ZtautauNp3 ZtautauNp4 ZtautauNp5

    # Z + bb + jets
    #ZeebbNp0 ZeebbNp1 ZeebbNp2 ZeebbNp3 ZeebbNp4 ZeebbNp5
    #ZmumubbNp0 ZmumubbNp1 ZmumubbNp2 ZmumubbNp3 ZmumubbNp4 ZmumubbNp5
    #ZtautaubbNp0 ZtautaubbNp1 ZtautaubbNp2 ZtautaubbNp3 ZtautaubbNp4 ZtautaubbNp5

    # Z + cc + jets
    #ZeeccNp0 ZeeccNp1 ZeeccNp2 ZeeccNp3 ZeeccNp4 ZeeccNp5
    #ZmumuccNp0 ZmumuccNp1 ZmumuccNp2 ZmumuccNp3 ZmumuccNp4 ZmumuccNp5
    #ZtautauccNp0 ZtautauccNp1 ZtautauccNp2 ZtautauccNp3 ZtautauccNp4 ZtautauccNp5

    # Core Samples
    # Low mass Z
    ZeeNp0Excl_Mll10to60 ZeeNp1Excl_Mll10to60 ZeeNp2Excl_Mll10to60 ZeeNp3Excl_Mll10to60 ZeeNp4Excl_Mll10to60
    ZmumuNp0Excl_Mll10to60 ZmumuNp1Excl_Mll10to60 ZmumuNp2Excl_Mll10to60 ZmumuNp3Excl_Mll10to60 ZmumuNp4Excl_Mll10to60
    ZtautauNp0Excl_Mll10to60 ZtautauNp1Excl_Mll10to60 ZtautauNp2Excl_Mll10to60 ZtautauNp3Excl_Mll10to60

    # W+jets
    #Sherpa_CT10_Wenu Sherpa_CT10_Wmunu Sherpa_CT10_Wtaunu

    # Core Sample
    # single top
    SingleTopWtChanIncl

    # ttbar
    #ttbar_LeptonFilter
    #PowhegPythia_AUET2BCT10_ttbar_LeptonFilter_AF2

    # Core Samples
    ttbarZ ttbarZj ttbarW ttbarWj
    TtbarHadTauhad
    TtbarLeptHad
    TtbarLeptLept
    TtbarLeptTauhad
    TtbarLeptTaulept
    TtbarTauhadTauhad
    TtbarTauleptHad
    TtbarTauleptTauhad
    TtbarTauleptTaulept

    # Core Samples
    # diboson
    llll_ZZ lllnu_WZ llnunu_ZZ llnunu_WW
    WpWmenuenu WpWmenumunu WpWmenutaunu
    WpWmmunuenu WpWmmunumunu WpWmmunutaunu
    WpWmtaunuenu WpWmtaunumunu WpWmtaunutaunu
    ZZ4lep ZZ4e ZZ4mu ZZ2e2mu
    llnunu_SS_EW6 llnunujj_SS

    enugammaPt10 munugammaPt10 taunugammaPt10
    eegammaPt10 mumugammaPt10 tautaugammaPt10
    ZWWStar_lllnulnu
    ZZZStar_nunullll
    WWWStar_lnulnulnu

    # HF
    #bbTomu20 ccTomu20

    # signal
    wA_noslep_WH_2Lep_2 
    wA_noslep_WH_2Lep_3 
    wA_noslep_WH_2Lep_5 
    wA_noslep_WH_2Lep_8 
    wA_noslep_WH_2Lep_9 
    wA_noslep_WH_2Lep_11
    wA_noslep_WH_2Lep_13
    wA_noslep_WH_2Lep_14
    wA_noslep_WH_2Lep_16
    wA_noslep_WH_2Lep_17
    wA_noslep_WH_2Lep_18
    wA_noslep_WH_2Lep_22
    wA_noslep_WH_2Lep_26
    wA_noslep_WH_2Lep_29
    wA_noslep_WH_2Lep_30
    wA_noslep_WH_2Lep_31
    wA_noslep_WH_2Lep_36
    wA_noslep_WH_2Lep_38
    wA_noslep_WH_2Lep_47
    wA_noslep_WH_2Lep_51
    wA_noslep_WH_2Lep_53
    wA_noslep_WH_2Lep_56
    wA_noslep_WH_2Lep_59
    wA_noslep_WH_2Lep_61
)

write_script () {
    local dset=$1
    local suffix=$2
    local options=""
    local template="batchSPSub.sh"
    local outFile="batchScripts/${dset}.sh"
    local inp=${dset}
    local out=${dset}${suffix}
    local opt=${options}
    cat  "${template}" \
        | sed "s/\${inp}/${inp}/g" \
        | sed "s/\${out}/${out}/g" \
        | sed "s/\${opt}/${opt}/g" > ${outFile}
    echo ${outFile}
    }

for d in ${DATA[@]}; do
    #echo "qsub -j oe -V -v inp=${d},out=${d}${append},opt=${option} -N ${d} -o batchlog batchSPSub.sh"
    scriptFile=$(write_script ${d} ${append})
    qsub -j oe -V -N ${d}_${append} -o batchLog ${scriptFile}
done

for d in ${MC[@]}; do
    #echo "qsub -j oe -V -v inp=${d},out=${d}${append},opt=${option} -N ${d}_${append} -o batchlog batchSPSub.sh"
    scriptFile=$(write_script ${d} ${append})
    qsub -j oe -V -N ${d}_${append} -o batchLog ${scriptFile}
done

echo "Finished submitting Susy Plot"
echo

