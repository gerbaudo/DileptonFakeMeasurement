#!/bin/bash


append="_AE_Dec17_n0115"
append="_Jan15_n0115"


# only one option supported: --1fb to use the 1/fb weights
#option="--1fb"
#option="--AD"

echo "Appending " ${append} " Run option: " ${option}

DATA=(

    # Data
    #periodA_Egamma periodA_Muons
    #periodB_Egamma periodB_Muons
    #periodC_Egamma periodC_Muons
    #periodD_Egamma periodD_Muons
    #periodE_Egamma periodE_Muons

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

)

for d in ${DATA[@]}; do

    qsub -j oe -V -v inp=${d},out=${d}${append},opt=${option} -N ${d} -o batchlog batchSPSub.sh
    #nohup SusyPlot -f filelist/${d}.txt -s ${d}${append} ${option} > batchLog/${d}${append}.susyplot.log &

done

for d in ${MC[@]}; do

    qsub -j oe -V -v inp=${d},out=${d}${append},opt=${option} -N ${d}_${append} -o batchlog batchSPSub.sh
    #qsub -j oe -V -v inp=${d},out=${d}${append},opt=${option} -N ${d}_${append} -o batchlog batchCFSub.sh

done

echo "Finished submitting Susy Plot"
echo

