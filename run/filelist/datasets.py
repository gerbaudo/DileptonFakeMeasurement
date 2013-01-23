#!/bin/env python

# List of datastets organized in a dict (imported from Matt, Jan 2013)
#
# davide.gerbaudo@gmail.com
# Jan 2013

# What files we want, with appropriate names
wantedDsets = { # mode : [dsets]
    'data' : []
    + ["period%(p)s.physics_%(s)s" % {'p':p, 's':s}
       for p in ['A', 'B', 'C', 'D', 'E'] for s in ['Egamma', 'Muons']]
    ,
    'mc12' : []
    ## # Alternative Z+jets AFII samples
    ## + ["Z%(ll)s%(f)%sJets_AF2" % {'ll':ll, 'f':f}
    ##    for ll in ['ee', 'mumu', 'tautau'] for f in ['Heavy', 'Light']]
    # New alternative Z+jets (Yippeee)
    + ["Sherpa_CT10_Z%s" % ll for ll in ['ee', 'mumu', 'tautau']]
    ##- # Z+jets
    ##- + ["Z%(ll)sNp%(np)d" % {'ll':ll, 'np':np}
    ##-    for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3, 4, 5]]
    ##- # Zbb + jets
    ##- + ["Z%(ll)sbbNp%(np)d" % {'ll':ll, 'np':np}
    ##-    for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3]]
    ##- # Zcc + jets
    ##- + ["Z%(ll)sccNp%(np)d" % {'ll':ll, 'np':np}
    ##-    for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3]]
    # Low mass Z
    + ["Z%(ll)sNp%(np)dExcl_Mll10to60" % {'ll':ll, 'np':np}
       for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3, 4]]
    ##- # W+Jets (temporary due to bugs)
    ##- + ["Sherpa_CT10_W%s" % lv for lv in ['enu', 'munu', 'taunu']]
    ##- # W+jets
    ##- + ["W%(lv)sNp%(np)d" % {'lv':lv, 'np':np}
    ##-    for lv in ['enu', 'munu', 'taunu'] for np in [0, 1, 2, 3, 4, 5]]
    ##- # Wbb
    ##- + ["WbbNp%d" % np for np in [0, 1, 2, 3]]
    ##- # Wcc
    ##- + ["WccNp%d" % np for np in [0, 1, 2, 3]]
    ##- + ["WcNp%d" % np for np in [0, 1, 2, 3]]
    # single top
    ##- + ["singletop_tchan_%s" % l for l in ['e', 'mu', 'tau']]
    + ["SingleTopWtChanIncl"]
    # ttbar
    ##- + ["McAtNloJimmy_CT10_ttbar_LeptonFilter",]
    + ["Ttbar%s" % ttd for ttd in ["LeptLept", "LeptTaulept", "TauleptTaulept",
                                   "LeptHad", "LeptTauhad", "TauleptHad",
                                   "TauleptTauhad", "HadTauhad", "TauhadTauhad",]]
    # + ["PowhegPythia_AUET2BCT10_ttbar_LeptonFilter_AF2",]
    + ["ttbar%s" % ttX for ttX in ["Z", "Zj", "W","Wj",]]
    # diboson
    + ["lllnu_WZ", "llll_ZZ", "llnunu_ZZ", "llnunu_WW",]
    + ["WpWm%s%s" % (lv1, lv2)
       for lv1 in ['enu', 'munu', 'taunu'] for lv2 in ['enu', 'munu', 'taunu']]
    + ["%sgammaPt10" % lv for lv in  ['enu', 'munu', 'taunu',]]
    + ["%sgammaPt10" % ll for ll in ['ee', 'mumu', 'tautau',]]
    # MISSING DIBOSON!!!!!
    + ["ZZ%s" % l4 for l4 in ['4lep', '4e', '4mu', '2e2mu',]]
    + ['llnunu_SS_EW6', 'llnunujj_SS',]
    # Triboson
    + ['ZWWStar_lllnulnu', 'ZZZStar_nunullll', 'WWWStar_lnulnulnu',]
    ##- # HF samples
    ##- + ["%(qq)sTo%(l)s15" % {'qq':qq, 'l':l} for qq in ['bb', 'cc'] for l in ['e', 'mu']]
    ,
    'susy' : []
    + ["wA_noslep_WH_2Lep_%d" % i for i in range(1, 61+1)]  # 2l
    + ["wA_noslep_WH_3Lep_%d" % i for i in range(1, 66+1)]  # 3l
    }
