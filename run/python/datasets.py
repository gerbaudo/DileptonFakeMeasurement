#!/bin/env python

# List of dataset definitions
#
# Details:
# see the Dataset class for more info on their attributes.
# Datasets can also be defined without being used; search for
# 'placeholder' to see the disabled samples.
#
# davide.gerbaudo@gmail.com
# Jan 2013

import unittest

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
    @property
    def isNotToBeMerged(self) : return self.isSignal
    @property
    def isToBeMerged(self) : return not self.isNotToBeMerged
    @property
    def isHeavyFlavor(self) : return self.group is 'heavyflavor'
    @property
    def isSignal(self) :
        return self.type is 'mc' and any(s in self.name for s in ('WH_2Lep', 'WH_2Lep'))
    @property
    def isSignalOrHiggs(self) : return self.isSignal or self.group is 'higgs'
    @property
    def isMcBackground(self) :
        "not written in stone (include higgs?), but that's what we need for the fake estimate"
        return self.type is 'mc' and not self.isSignalOrHiggs

def allGroups(datasets=[]) : return list(set(d.group for d in datasets))
def allDatasets(datasets=[]) : return list(set(d.name for d in datasets))
def activeDatasets(datasets=[]) : return filter(lambda d : not d.placeholder, datasets)
def setSameGroupForAllData(datasets=[], group='data') :
    for d in datasets :
        if d.type is 'data' : d.group = group
    return datasets

datasets = []
sampleType, group, process = None, None, None
placeholder = True


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

group = 'zjets'
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
nps = ['0Excl', '1Excl', '2Excl', '3Excl', '4Excl', '5Incl']
template, process = "AlpgenJimmy_Auto_AUET2CTEQ6L1_Z%(ll)sNp%(np)s_Mll10to60", 'Zlljets_Mll10to60'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll, 'np':np}, process)
             for ll, dsids in [('ee', range(146830, 146835+1))]
             for d, np in rzip(dsids, nps)]
nps = ['0Excl', '1Excl', '2Excl', '3Excl', '4Excl', '5Incl']
template, process = "AlpgenJimmy_Auto_Z%(ll)sNp%(np)s_Mll10to60", 'Zlljets_Mll10to60'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll, 'np':np}, process)
             for ll, dsids in [('mumu', range(146840,   146845+1)),
                               ('tautau', range(146850, 146855+1))]
             for d, np in rzip(dsids, nps)]

group = 'wjets'
template, process = "Sherpa_CT10_W%(lv)s", 'Sherpa W+jets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv}, process, placeholder)
             for lv, d in [('enu',   147774),
                           ('munu',  147775),
                           ('taunu', 147776)]]
nps = [0, 1, 2, 3, 4, 5]
template, process = "AlpgenJimmy_AUET2CTEQ6L1_W%(lv)sNp%(np)d", 'Wlvjets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv, 'np':np}, process,
                     placeholder)
             for lv, dsids in [('enu',   range(107680, 107685+1)),
                               ('munu',  range(107690, 107695+1)),
                               ('taunu', range(107700, 107705+1))]
             for d, np in rzip(dsids, nps)]
nps = [0, 1, 2, 3, 4, 5]
template, process = "AlpgenPythia_P2011C_W%(lv)sNp%(np)d", 'Wlvjets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv, 'np':np}, process)
             for lv, dsids in [('enu',   range(117680, 117685+1)),
                               ('munu',  range(117690, 117695+1)),
                               ('taunu', range(117700, 117705+1))]
             for d, np in rzip(dsids, nps)]
nps = ['0', '1', '2', '3', '4', '5incl']
template, process = "AlpgenPythia_Auto_P2011C_W%(lv)sNp%(np)s", 'Wlvjets'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv, 'np':np}, process,
                     placeholder)
             for lv, dsids in [('enu',   range(147025, 147030+1)),
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
datasets += [Dataset(sampleType, d, group, template%{'qq':qq, 'np':np}, process)
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
template, process = 'McAtNloJimmy_AUET2CT10_SingleTopSChanW%(lv)s', 'singletop'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv}, process)
             for d, lv in [(108343, 'enu'), (108344, 'munu'), (108345, 'taunu')]]
name, process = 'McAtNloJimmy_AUET2CT10_SingleTopWtChanIncl', 'singletop'
datasets += [Dataset(sampleType, 108346, group, name, process)]
template, process = 'AcerMCPythia_AUET2BCTEQ6L1_singletop_tchan_%(l)s', 'singletop'
datasets += [Dataset(sampleType, d, group, template%{'l':l}, process)
             for d, l in [(117360, 'e'), (117361, 'mu'), (117362, 'tau')]]
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
template, process = "PowhegPythia8_AU2CT10_WpWm_%(ll)s", 'PowhegPythia8_WpWm'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll}, process)
             for d, ll in [(126928, 'ee'), (126929, 'me'), (126930, 'te'),
                           (126931, 'em'), (126932, 'mm'), (126933, 'tm'),
                           (126934, 'et'), (126935, 'mt'), (126936, 'tt')]]
template, process = "gg2wwJimmy_AUET2CT10_WpWm%(lvlv)s", 'gg2wwJimmy'
datasets += [Dataset(sampleType, d, group, template%{'lvlv':lvlv}, process)
             for d, lvlv in [(169471, 'enuenu'),     (169472, 'enumunu'),  (169473, 'enutaunu'),
                             (169474, 'munumunu'),   (169475, 'munuenu'),  (169476, 'munutaunu'),
                             (169477, 'taunutaunu'), (169478, 'taunuenu'), (169479, 'taunumunu')]]
template, process = "gg2ZZJimmy_AUET2CT10_ZZ%(l4)s", 'gg2ZZJimmy'
datasets += [Dataset(sampleType, d, group, template%{'l4':l4}, process)
             for d, l4 in [(116600, '4lep'),
                           (116601, '4e'),
                           (116602, '4mu'),
                           (116603, '2e2mu')]]
template, process = "Sherpa_CT10_%(lepVV)s", 'Sherpa_VV_lep'
datasets += [Dataset(sampleType, d, group, template%{'lepVV':lepVV}, process, placeholder)
             for d, lepVV in [(126892, 'llnunu_WW'), (126893, 'lllnu_WZ'),
                              (161963, 'llnunu_ZZ'), (174834, 'llll_ZZ')]]
template, process = "Sherpa_CT10_%(lv)sgammaPt10", 'Sherpa_Vgamma'
datasets += [Dataset(sampleType, d, group, template%{'lv':lv}, process)
             for d, lv in [(126739, 'enu'),   (126742, 'munu'), (126856, 'taunu')]]
template, process = "Sherpa_CT10_%(ll)sgammaPt10", 'Sherpa_Vgamma'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll}, process, placeholder)
             for d, ll in [(145161, 'ee'),    (145162, 'mumu'), (126854, 'tautau')]]
template, process = "Sherpa_CT10_%(llss)s", 'Sherpa_llnunu'
datasets += [Dataset(sampleType, d, group, template%{'llss':llss}, process)
             for d, llss in [(126988, 'llnunu_SS_EW6'), (126989, 'llnunujj_SS')]]
template, process = "Sherpa_CT10_VVto%(l)snuqq", 'Sherpa_VVtolnuqq'
datasets += [Dataset(sampleType, d, group, template%{'l':l}, process)
             for d, l in [(157817, 'e'), (157818, 'mu'), (157819, 'tau')]]
template, process = "Sherpa_CT10_VVto%(ll)sqq", 'Sherpa_VVtollqq'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll}, process)
             for d, ll in [(157814, 'ee'), (157815, 'mumu'), (157816, 'tautau')]]
template, process = "MadGraphPythia_AUET2BCTEQ6L1_%(VVV)sStar_%(fs)s", 'Triboson'
datasets += [Dataset(sampleType, d, group, template%{'VVV':VVV, 'fs':fs}, process)
             for d, VVV, fs in [(167006, 'WWW', 'lnulnulnu'),
                                (167007, 'ZWW', 'lllnulnu'),
                                (167008, 'ZZZ', 'nunullll')]]
template, process = "PowhegPythia8_AU2CT10_ZZ_%(l4)s_mll4_2pt5", 'PowhegPythia8_ZZ'
datasets += [Dataset(sampleType, d, group, template%{'l4':l4}, process)
             for d, l4 in [(126937, '4e'),  (126938, '2e2mu'),   (126939, '2e2tau'),
                           (126940, '4mu'), (126941, '2mu2tau'), (126942, '4tau'),]]
template, process = "PowhegPythia8_AU2CT10_ZZllnunu_%(ll)s_mll4", 'PowhegPythia8_ZZ'
datasets += [Dataset(sampleType, d, group, template%{'ll':ll}, process)
             for d, ll in [(126949, 'ee'), (126950, 'mm'), (126951, 'tt')]]
template, process = "PowhegPythia8_AU2CT10_WZ_%(lvll)s_mll%(mll)s_2L5", 'PowhegPythia8_WZ'
datasets += [Dataset(sampleType, d, group, template%{'lvll':lvll, 'mll':mll}, process)
             for d, lvll, mll in [(129477, 'Wm11Z11', '0p250d0'),
                                  (129478, 'Wm11Z13', '0p4614d0'),
                                  (129479, 'Wm11Z15', '3p804d0'),
                                  (129480, 'Wm13Z11', '0p250d0'),
                                  (129481, 'Wm13Z13', '0p4614d0'),
                                  (129482, 'Wm13Z15', '3p804d0'),
                                  (129483, 'Wm15Z11', '0p250d0'),
                                  (129484, 'Wm15Z13', '0p4614d0'),
                                  (129485, 'Wm15Z15', '3p804d0'),
                                  (129486, 'W11Z11',  '0p250d0'),# these are actually Wp
                                  (129487, 'W11Z13',  '0p4614d0'),
                                  (129488, 'W11Z15',  '3p804d0'),
                                  (129489, 'W13Z11',  '0p250d0'),
                                  (129490, 'W13Z13',  '0p4614d0'),
                                  (129491, 'W13Z15',  '3p804d0'),
                                  (129492, 'W15Z11',  '0p250d0'),
                                  (129493, 'W15Z13',  '0p4614d0'),
                                  (129494, 'W15Z15',  '3p804d0')]]

group = 'heavyflavor'
template, process = "Pythia8B_AU2_CTEQ6L1_%(qq)sTomu20", 'HF qqbar'
datasets += [Dataset(sampleType, d, group, template%{'qq':qq}, process)
             for d, qq in [(129136, 'bb'), (147668, 'cc')]]

group = 'higgs'
template, process = "PowhegPythia8_AU2CT10_ggH125_%(fs)s", 'ggH125'
datasets += [Dataset(sampleType, d, group, template%{'fs':fs}, process)
             for d, fs in [(160155, 'ZZ4lep'), (160655, 'ZZllnunu'),
                           (161005, 'WW2lep_EF_15_5')]]
template, process = "PowhegPythia8_AU2CT10_VBFH125_%(fs)s", 'VBFH125'
datasets += [Dataset(sampleType, d, group, template%{'fs':fs}, process)
             for d, fs in [(160205, 'ZZ4lep'), (160705, 'ZZllnunu'),
                           (161055, 'WW2lep_EF_15_5')]]
template, process = "Pythia8_AU2CTEQ6L1_WH125_%(fs)s", 'WH125'
datasets += [Dataset(sampleType, d, group, template%{'fs':fs}, process)
             for d, fs in [(160255, 'ZZ4lep'), (160505, 'ZZllqq'), (160755, 'ZZllnunu'),
                           (161105, 'WW2lep'), (161805, 'lnubb')]]
template, process = "Pythia8_AU2CTEQ6L1_ZH125_%(fs)s", 'ZH125'
datasets += [Dataset(sampleType, d, group, template%{'fs':fs}, process)
             for d, fs in [(160305, 'ZZ4lep'), (160555, 'ZZllqq'), (160805, 'ZZllnunu'),
                           (161155, 'WW2lep'),
                           (161675, 'tautaull'), (161686, 'tautaulh'), (161697, 'tautauhh'),
                           (167418, 'mumu')]]
template, process = "Pythia8_AU2CTEQ6L1_ttH125_%(fs)s", 'ttH125'
datasets += [Dataset(sampleType, d, group, template%{'fs':fs}, process)
             for d, fs in [(161305, 'WWinclusive'),
                           (161708, 'tautaull'), (161719, 'tautaulh'), (161730, 'tautauhh'),
                           (169072, 'ZZinclusive')]]

groupTemplate = 'WH_2Lep_%d'
template, process = "Herwigpp_simplifiedModel_wA_noslep_WH_2Lep_%d", "wA_noslep_WH_2Lep_%d"
datasets += [Dataset(sampleType, d, groupTemplate%nth, template%nth, process%nth)
             for d, nth in rzip(range(176574, 176634+1), range(1, 61+1))]
groupTemplate = 'WH_3Lep_%d'
template, process = "Herwigpp_simplifiedModel_wA_noslep_WH_3Lep_%d", "wA_noslep_WH_3Lep_%d"
datasets += [Dataset(sampleType, d, groupTemplate%nth, template%nth, process%nth, placeholder)
             for d, nth in rzip(range(176641, 176706+1), range(1, 66+1))]
groupTemplate = 'notauhad_WH_2Lep_%d'
template, process = "Herwigpp_sM_wA_noslep_notauhad_WH_2Lep_%d", "wA_noslep_notauhad_WH_2Lep_%d"
datasets += [Dataset(sampleType, d, groupTemplate%nth, template%nth, process%nth)
             for d, nth in rzip(range(177501, 177528+1), range(1, 28+1))]


#
# testing
#

class CategorizationTest(unittest.TestCase) :
    def testSignalIsConsistent(self) :
        for d in datasets :
            notData = d.type is not 'data'
            notBkg = d.type is 'mc' and d.isSignal
            notHf = not d.isHeavyFlavor
            isSignal = d.isSignal
            self.assertEqual(isSignal, (notData and notBkg and notHf),
                             d.name+' : '
                             +', '.join(["%s : %s"%(k, eval(k))
                                         for k in ['isSignal', 'notData', 'notBkg', 'notHf']]))

#
# testing
#

class CategorizationTest(unittest.TestCase) :
    def testSignalIsConsistent(self) :
        for d in datasets :
            notData = d.type is not 'data'
            notBkg = d.type is 'mc' and d.isSignal
            notHf = not d.isHeavyFlavor
            isSignal = d.isSignal
            self.assertEqual(isSignal, (notData and notBkg and notHf),
                             d.name+' : '
                             +', '.join(["%s : %s"%(k, eval(k))
                                         for k in ['isSignal', 'notData', 'notBkg', 'notHf']]))

if __name__=='__main__' :
    def filterByGroup(dsets) :
        groups = sorted(d.group for d in dsets)
        return dict((g, filter(lambda x: g==x.group, dsets)) for g in groups)
    def filterByProcess(dsets) :
        processes = sorted(d.process for d in dsets)
        return dict((p, filter(lambda x: p==x.process, dsets)) for p in processes)
    def printByGroupByProcess(dsets) :
        for g, gdsets in filterByGroup(dsets).iteritems() :
            print "------- %s ------"%g
            for p, pdsets in filterByProcess(gdsets).iteritems() :
                print "---   %s    ---"%p
                print "Group Dsid Name Type Process"
                print '\n'.join(["%s %s %s %s %s"%(s.group, s.dsid, s.name, s.type,
                                                   s.process)
                                 for s in pdsets])
    unusedDsets = filter(lambda x:     x.placeholder, datasets)
    usedDsets   = filter(lambda x: not x.placeholder, datasets)
    linebreak = '-'*20
    print "List of available datasets[%d]"%len(datasets)
    print linebreak
    print "Unused datasets (placeholders)"
    print linebreak
    printByGroupByProcess(unusedDsets)
    print linebreak
    print linebreak
    print "Used datasets"
    print linebreak
    printByGroupByProcess(usedDsets)
    print linebreak
    print
    print
    print "Summary:"
    print "Total number of datasets : %d"%len(datasets)
    print "  unused : %d"%len(unusedDsets)
    print '\n'.join("         : %d : %s"%(len(gdsets), g)
                    for g, gdsets in filterByGroup(unusedDsets).iteritems())
    print "    used : %d"%len(usedDsets)
    print '\n'.join("         : %d : %s"%(len(gdsets), g)
                    for g, gdsets in filterByGroup(usedDsets).iteritems())
    print linebreak
    print 'Testing...'
    unittest.main()
