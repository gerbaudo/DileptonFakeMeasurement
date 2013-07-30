#!/bin/env python

# List of datastets organized in a dict (imported from Matt, Jan 2013)
#
# davide.gerbaudo@gmail.com
# Jan 2013

# What files we want, with appropriate names
wantedDsets = { # mode : [dsets]
    'data' : []
    + ["period%(p)s.physics_%(s)s" % {'p':p, 's':s}
       for p in ['A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'L'] for s in ['Egamma', 'Muons']]
    ,
    'mc12' : []
    ## # Alternative Z+jets AFII samples
    ## + ["Z%(ll)s%(f)%sJets_AF2" % {'ll':ll, 'f':f}
    ##    for ll in ['ee', 'mumu', 'tautau'] for f in ['Heavy', 'Light']]
    ## # New alternative Z+jets (Yippeee)
    ## + ["Sherpa_CT10_Z%s" % ll for ll in ['ee', 'mumu', 'tautau']]
    # Z+jets
    + ["AlpgenPythia_P2011C_Z%(ll)sNp%(np)d" % {'ll':ll, 'np':np}
       for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3, 4, 5]]
    # Zbb + jets
    + ["AlpgenPythia_P2011C_Z%(ll)sbbNp%(np)d" % {'ll':ll, 'np':np}
       for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3]]
    # Zcc + jets
    + ["AlpgenPythia_P2011C_Z%(ll)sccNp%(np)d" % {'ll':ll, 'np':np}
       for ll in ['ee', 'mumu', 'tautau'] for np in [0, 1, 2, 3]]
    # Low mass Z
    + ["AlpgenPythia_P2011C_Z%(ll)sNp%(np)dExcl" % {'ll':ll, 'np':np}
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
    + ['McAtNloJimmy_CT10_ttbar_LeptonFilter']
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

def rzip(*iterables) :
    """rigid zip: ensures input iterables have the same length.
    Used to make sure that the dsid ranges match the dset names."""
    lengths = list(set([len(l) for l in iterables]))
    assert len(lengths)==1,"zip lists with etherogeneus lengths \n"\
           "%s"%'\n'.join(["[%d] : %s"%(len(l), str(l)) for l in iterables])
    return zip(*iterables)

class Dataset :
    """Container to uniquely specify a dataset (through a dsid), and
    specify additional user-friendly attributes"""
    def __init__(self, sampleType, dsid=None, group=None, name=None, process=None, placeholder=False) :
        assert sampleType in ['data','mc'], "sampleType: %s for %s"%(sampleType, name)
        self.type = sampleType # data or mc
        self.dsid = dsid # 6-digit id for mc; none for data
        self.group = group # group in which the histo will appear in the stack
        self.name = name # full name, usually as in <stuff>.<dsid>.<name>.SusyNt.<tags>
        self.process = process # physical process, for example 'Zbb + jets' (short, generic, usually appears as the common root of the name)
        self.placeholder = placeholder # just a placeholder, its job won't be submitted


datasets = []
sampleType, group, process = None, None, None
placeholder = True # search for 'placeholder' to see the disabled samples


sampleType = 'data'
periods = ['A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'L']
group, dsid = 'Egamma', None
datasets += [Dataset(sampleType, dsid, group, n, process)
             for n in ["period%(p)s.physics_%(g)s.PhysCont" % {'p':p, 'g':group}
                       for p in periods]]
group, dsid = 'Muons', None
datasets += [Dataset(sampleType, dsid, group, n, process)
             for n in ["period%(p)s.physics_%(g)s.PhysCont" % {'p':p, 'g':group}
                       for p in periods]]

sampleType = 'mc'

group = 'Zjets'
nps = [0, 1, 2, 3, 4, 5]
template, process = "AlpgenPythia_P2011C_Z%(ll)sNp%(np)d", 'Zlljets'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll, 'np':np}, process)
             for ll, dsids in [('ee',     range(117650, 117655+1)),
                               ('mumu',   range(117660, 117665+1)),
                               ('tautau', range(117670, 117675+1))]
             for d, np in rzip(dsids, nps)]
nps = [0, 1, 2, 3]
template, process = "AlpgenPythia_P2011C_Z%(ll)sbbNp%(np)d", 'Zbbjets'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll, 'np':np}, process)
             for ll, dsids in [('ee',     range(110817, 110820+1)),
                               ('mumu',   range(110821, 110824+1)),
                               ('tautau', range(110825, 110828+1))]
             for d, np in rzip(dsids, nps)]
nps = [0, 1, 2, 3]
template, process = "AlpgenPythia_P2011C_Z%(ll)sccNp%(np)d", 'Zccjets'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll, 'np':np}, process)
             for ll, dsids in [('ee',     range(110805, 110808+1)),
                               ('mumu',   range(110809, 110812+1)),
                               ('tautau', range(110813, 110816+1))]
             for d, np in rzip(dsids, nps)]
template, process = "Sherpa_CT10_Z%(ll)s", 'Sherpa Z+jets'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll}, process, placeholder)
             for ll, d in [('ee',     147770),
                           ('mumu',   147771),
                           ('tautau', 147772)]]

group = 'Wjets'
template, process = "Sherpa_CT10_W%(lv)s", 'Sherpa W+jets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv}, process, placeholder)
             for ll, d in [('enu',   147774),
                           ('munu',  147775),
                           ('taunu', 147776)]]
nps = [0, 1, 2, 3, 4, 5]
template, process = "AlpgenJimmy_AUET2CTEQ6L1_W%(lv)sNp%(np)d", 'Wlvjets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv, 'np':np}, process,
                     placeholder)
             for ll, dsids in [('enu',   range(107680, 107685+1)),
                               ('munu',  range(107690, 107695+1)),
                               ('taunu', range(107700, 107705+1))]
             for d, np in rzip(dsids, nps)]
nps = [0, 1, 2, 3, 4, 5]
template, process = "AlpgenPythia_P2011C_W%(lv)sNp%(np)d", 'Wlvjets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv, 'np':np}, process,
                     placeholder)
             for ll, dsids in [('enu',   range(117680, 117685+1)),
                               ('munu',  range(117690, 117695+1)),
                               ('taunu', range(117700, 117705+1))]
             for d, np in rzip(dsids, nps)]
nps = ['0', '1', '2', '3', '4', '5incl']
template, process = "AlpgenPythia_Auto_P2011C_W%(lv)sNp%(np)s", 'Wlvjets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv, 'np':np}, process,
                     placeholder)
             for ll, dsids in [('enu',   range(147025, 147030+1)),
                               ('munu',  range(147033, 147038+1)),
                               ('taunu', range(147041, 147046+1))]
             for d, np in rzip(dsids, nps)]
npsBase = [0, 1, 2, 3]
template, process = "AlpgenJimmy_AUET2CTEQ6L1_W%(qq)sNp%(np)d", 'Wqqjets'
datasets += [Dataset(sampleType, d, group, template%{'qq':qq, 'np':np}, process,
                     placeholder)
             for qq, dsids, nps in [('bb', range(107280, 107283+1), npsBase),
                                    ('cc', range(117284, 117287+1), npsBase),
                                    ('c',  range(117293, 117297+1), npsBase+[4])]
             for d, np in rzip(dsids, nps)]
npsBase = [0, 1, 2, 3]
template, process = "AlpgenPythia_P2011C_W%(qq)sNp%(np)d", 'Wqqjets'
datasets += [Dataset(sampleType, d, group, template%{'qq':qq, 'np':np}, process,
                     placeholder)
             for qq, dsids, nps in [('bb', range(110801, 110804+1), npsBase),
                                    ('cc', range(126606, 126609+1), npsBase),
                                    ('c',  range(126601, 126605+1), npsBase+[4])]
             for d, np in rzip(dsids, nps)]
npsBase = ['0', '1', '2', '3incl']
template, process = "AlpgenPythia_Auto_P2011C_W%(qq)sNp%(np)s", 'Wqqjets'
datasets += [Dataset(sampleType, d, group, template%{'qq':qq, 'np':np}, process,
                     placeholder)
             for qq, dsids, nps in [('bb', range(200256, 200259+1), npsBase),
                                    ('cc', range(200156, 200159+1), npsBase),
                                    ('c',  range(200056, 200060+1),
                                     npsBase[:-1]+['3', '4incl'])]
             for d, np in rzip(dsids, nps)]

group = 'ttbar'
name, process = 'McAtNloJimmy_AUET2CT10_SingleTopWtChanIncl', 'singletop'
datasets += [Dataset(sampleType, 108346, group, name, process)]
name, process = 'McAtNloJimmy_CT10_ttbar_LeptonFilter', 'ttbar'
datasets += [Dataset(sampleType, 105200, group, name, process)]
template, process = "MadGraphPythia_AUET2BCTEQ6L1_ttbar%(ttX)s", 'ttbarV'
datasets += [Dataset(sampleType, d, group, template%{'ttX':ttX}, process)
             for d, ttX in [(119353, 'W' ), (119354, 'Wj'),
                            (119355, 'Z' ), (119356, 'Zj')]]
name, process = 'MadgraphPythia_AUET2B_CTEQ6L1_ttbarWW', 'ttbarWW'
datasets += [Dataset(sampleType, 119583, group, name, process)]

#- ["singletop_tchan_%s" % l for l in ['e', 'mu', 'tau']]
#- ["Ttbar%s" % ttd for ttd in ["LeptLept", "LeptTaulept", "TauleptTaulept",
#-                              "LeptHad", "LeptTauhad", "TauleptHad",
#-                              "TauleptTauhad", "HadTauhad", "TauhadTauhad",]]
#- ["PowhegPythia_AUET2BCT10_ttbar_LeptonFilter_AF2",]

group = 'diboson'

template, process = "gg2wwJimmy_AUET2CT10_WpWm%(lvlv)s", 'gg2wwJimmy'
datasets += [Dataset(sampleType, d, group, template%{'lvlv':lvlv}, process)
             for d, lvlv in [(169471, 'WpWmenuenu'), (169472, 'WpWmenumunu'), (169473, 'WpWmenutaunu'),
                             (169474, 'WpWmmunumunu'), (169475, 'WpWmmunuenu'), (169476, 'WpWmmunutaunu'),
                             (169477, 'WpWmtaunutaunu'), (169478, 'WpWmtaunuenu'), (169479, 'WpWmtaunumunu')]]
template, process = "gg2ZZJimmy_AUET2CT10_ZZ%(l4)s", 'gg2ZZJimmy'
datasets += [Dataset(sampleType, d, group, template%{'l4':l4}, process)
             for d, l4 in [(116601, '4e'),
                           (116602, '4mu'),
                           (116603, '2e2mu')]]
template, process = "Sherpa_CT10_%(lepVV)s", 'Sherpa_VV_lep'
datasets += [Dataset(sampleType, d, group, template%{'lepVV':lepVV}, process)
             for d, lepVV in [(174834, 'llll_ZZ'),   (161963, 'llnunu_ZZ'),
                              (126892, 'llnunu_WW'), (161961, 'lllnu_WZ')]]
template, process = "Sherpa_CT10_%(ll)sPt10", 'Sherpa_Vgamma'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll}, process)
             for d, ll in [(145161, 'eegamma'),   (145162, 'mumugamma'),
                           (126739, 'enugamma'),   (126742, 'munugamma'),
                           (126856, 'taunugamma'), (126854, 'tautaugamma')]]
template, process = "Sherpa_CT10_%(llss)s", 'Sherpa_llnunu'
datasets += [Dataset(sampleType, d, group, template%{'llss':llss}, process)
             for d, llss in [(126988, 'llnunu_SS_EW6'), (126989, 'llnunujj_SS')]]

#- 160305.Pythia8_AU2CTEQ6L1_ZH125_ZZ4lep
#- # Triboson
#- + ['ZWWStar_lllnulnu', 'ZZZStar_nunullll', 'WWWStar_lnulnulnu',]
#- # HF samples
#- # + ["%(qq)sTo%(l)s15" % {'qq':qq, 'l':l} for qq in ['bb', 'cc'] for l in ['e', 'mu']]

group = None
template, process = "Herwigpp_simplifiedModel_wA_noslep_WH_2Lep_%d", "wA_noslep_WH_2Lep_%d"
datasets += [Dataset(sampleType, d, group, template%nth, process%nth)
             for d, nth in rzip(range(176574, 176634+1), range(1, 61+1))]
template, process = "Herwigpp_simplifiedModel_wA_noslep_WH_3Lep_%d", "wA_noslep_WH_3Lep_%d"
datasets += [Dataset(sampleType, d, group, template%nth, process%nth)
             for d, nth in rzip(range(176641, 176706+1), range(1, 66+1))]


if __name__=='__main__' :
    print "List of available datasets[%d]"%len(datasets)
    print '-'*20
    print "Group Dsid Name Type Process"
    print '\n'.join(["%s %s %s %s %s"%(s.group, s.dsid, s.name, s.type, s.process) for s in datasets])
