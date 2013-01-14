# Make lists of input root files

import subprocess

# Directory where files are
#dir = "/gdata/atlas/ucintprod/SusyNt/mc12_n0041/"
#dir = "/gdata/atlas/ucintprod/SusyNt/mc12_n0105/"
#dir = "/gdata/atlas/ucintprod/SusyNt/mc12_n0111/"
#dir = "/gdata/atlas/ucintprod/SusyNt/mc12_n0114/"
dir = "/gdata/atlas/ucintprod/SusyNt/mc12_n0115/"

# Tag for production
tags = []
tags.append("n0115")
#tags.append("n0111")

# What files we want, with appropriate names
wanted = (

    # Alternative Z+jets AFII samples
    #"ZeeLightJets_AF2", "ZmumuLightJets_AF2", "ZtautauLightJets_AF2",
    #"ZeeHeavyJets_AF2", "ZmumuHeavyJets_AF2", "ZtautauHeavyJets_AF2",

    # New alternative Z+jets (Yippeee)
    "Sherpa_CT10_Zee",
    "Sherpa_CT10_Zmumu",
    "Sherpa_CT10_Ztautau",

    # Z+jets
    "ZeeNp0", "ZeeNp1", "ZeeNp2", "ZeeNp3", "ZeeNp4", "ZeeNp5",
    "ZmumuNp0", "ZmumuNp1", "ZmumuNp2", "ZmumuNp3", "ZmumuNp4", "ZmumuNp5",
    "ZtautauNp0", "ZtautauNp1", "ZtautauNp2", "ZtautauNp3", "ZtautauNp4", "ZtautauNp5",

    # Zbb + jets
    "ZeebbNp0", "ZeebbNp1", "ZeebbNp2", "ZeebbNp3",
    "ZmumubbNp0", "ZmumubbNp1", "ZmumubbNp2", "ZmumubbNp3",
    "ZtautaubbNp0", "ZtautaubbNp1", "ZtautaubbNp2", "ZtautaubbNp3",

    # Zcc + jets
    "ZeeccNp0", "ZeeccNp1", "ZeeccNp2", "ZeeccNp3",
    "ZmumuccNp0", "ZmumuccNp1", "ZmumuccNp2", "ZmumuccNp3",
    "ZtautauccNp0", "ZtautauccNp1", "ZtautauccNp2", "ZtautauccNp3",
    
    # Low mass Z
    "ZeeNp0Excl_Mll10to60", "ZeeNp1Excl_Mll10to60", "ZeeNp2Excl_Mll10to60",
    "ZeeNp3Excl_Mll10to60", "ZeeNp4Excl_Mll10to60",
    "ZmumuNp0Excl_Mll10to60", "ZmumuNp1Excl_Mll10to60", "ZmumuNp2Excl_Mll10to60",
    "ZmumuNp3Excl_Mll10to60", "ZmumuNp4Excl_Mll10to60",
    "ZtautauNp0Excl_Mll10to60", "ZtautauNp1Excl_Mll10to60", "ZtautauNp2Excl_Mll10to60",
    "ZtautauNp3Excl_Mll10to60", "ZtautauNp4Excl_Mll10to60", 

    # W+Jets (temporary due to bugs)
    "Sherpa_CT10_Wenu", "Sherpa_CT10_Wmunu", "Sherpa_CT10_Wtaunu",

    # W+jets
    "WenuNp0", "WenuNp1", "WenuNp2", "WenuNp3", "WenuNp4", "WenuNp5",
    "WmunuNp0", "WmunuNp1", "WmunuNp2", "WmunuNp3", "WmunuNp4", "WmunuNp5",
    "WtaunuNp0", "WtaunuNp1", "WtaunuNp2", "WtaunuNp3", "WtaunuNp4", "WtaunuNp5",

    # Wbb
    "WbbNp0", "WbbNp1", "WbbNp2", "WbbNp3",

    # Wcc
    "WccNp0", "WccNp1", "WccNp2", "WccNp3",
    "WcNp0", "WcNp1", "WcNp2", "WcNp3",

    # single top
    "singletop_tchan_e", "singletop_tchan_mu", "singletop_tchan_tau", "SingleTopWtChanIncl",

    # ttbar
    "McAtNloJimmy_CT10_ttbar_LeptonFilter",
    "TtbarLeptLept", "TtbarLeptTaulept", "TtbarTauleptTaulept",
    "TtbarLeptHad", "TtbarLeptTauhad", "TtbarTauleptHad",
    "TtbarTauleptTauhad", "TtbarHadTauhad", "TtbarTauhadTauhad",
    #"PowhegPythia_AUET2BCT10_ttbar_LeptonFilter_AF2",
    "ttbarZ", "ttbarZj",
    "ttbarW","ttbarWj",

    # diboson
    "lllnu_WZ", "llll_ZZ", "llnunu_ZZ", "llnunu_WW",
    "WpWmenuenu", "WpWmenumunu", "WpWmenutaunu",
    "WpWmmunuenu", "WpWmmunumunu", "WpWmmunutaunu",
    "WpWmtaunuenu", "WpWmtaunumunu", "WpWmtaunutaunu",
    "enugammaPt10", "munugammaPt10", "taunugammaPt10",
    "eegammaPt10", "mumugammaPt10", "tautaugammaPt10",

    # MISSING DIBOSON!!!!!
    "ZZ4lep", "ZZ4e", "ZZ4mu", "ZZ2e2mu",
    "llnunu_SS_EW6", "llnunujj_SS",

    # Triboson
    "ZWWStar_lllnulnu",
    "ZZZStar_nunullll",
    "WWWStar_lnulnulnu",

    # HF samples
    "bbTomu15", "bbToe15",
    "ccTomu15", "ccToe15",
    )


###############################################################################################
#                           Don't need to edit below here!!!                                  #
###############################################################################################
dlist = []
for tag in tags:
    print tag
    ls = subprocess.Popen(["ls " + dir + " | grep " + tag + " | grep user"],shell=True,stdout=subprocess.PIPE)
    templist = (ls.stdout.read()).split("\n")
    del templist[-1]
    for dset in templist:
        dlist.append( dset )
    

def contains(dataset, name):
    if (name + ".") in dataset:
        return True
    return False

def makeFile(dataset, name):
    ls = subprocess.Popen(["ls " + dir + dataset + "/* > " + name + ".txt"],shell=True)
    ls.wait()

for ds in dlist:
    print ds
    for name in wanted:
        if(contains(ds,name) and not ("_a" in ds)):
            makeFile(ds,name)

    
