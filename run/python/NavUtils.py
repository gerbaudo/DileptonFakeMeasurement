#!/bin/env python

# Utility classes and function to navigate the histograms in a root file
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, re, unittest
import ROOT as r

class HistoType(object):
    "Type of histogram, defined by plot region, channel, variable, syst"
    attributes = ['pr', 'ch', 'var', 'syst']
    def __init__(self, pr='', ch='', var='', syst=''):
        for att in HistoType.attributes : setattr(self, att, eval(att))
    def sameas(self, rhs):
        return all([getattr(self,att)==getattr(rhs,att) for att in HistoType.attributes])
    def __eq__(self, other) : return self.sameas(other)
    def __str__(self) : return ', '.join(["%s : %s"%(a, getattr(self,a)) for a in HistoType.attributes])
    def __hash__(self) : return hash(self.__str__())
    def matchAllAvailabeAttrs(self, rhs) :
        return all([ar==al for ar,al in [(getattr(self,a), getattr(rhs,a))
                                         for a in HistoType.attributes] if ar])

def setHistoType(h, type) : setattr(h, 'type', type)
def setHistoSample(h, sample) : setattr(h, 'sample', sample)

class HistoNameClassifier :
    """Extract pieces of information from the histogram name and attach an HistoType attribute
    See the NEWHIST macro in SusyPlotter::Begin
    """
    def __init__(self, verbose=False) :
        self.verbose = verbose
        self.rep = re.compile('(?P<pr>.*?)_'   # plot region (non greedy)
                              +'(?P<ch>.*?)_'  # channel (non greedy)
                              +'(?P<var>.*)_'  # var name (greedy, can contain '_')
                              +'(?P<syst>.*)') # last token, everything that's left
    def histoType(self, histoname='') :
        match = self.rep.search(histoname)
        if not match :
            if self.verbose : print "cannot classify %s" % histoname
        else :
            kargs = dict([(g, match.group(g)) for g in ['pr', 'ch', 'var', 'syst']])
            return  HistoType(**kargs)

def getAllHistoNames(inputDir, verbose=False, onlyTH1=False, onlyTH2=False, onlyTH3=False) :
    "Provide a list of all histograms in the file; search recursively (use FindObjectAny to retrieve from subdirs"
    assert onlyTH1 + onlyTH2 + onlyTH2 <= 1, "specify onlyTH* one at the time : %s" % ' '.join(["%s=%s"%(o, eval(o)) for o in ['onlyTH1', 'onlyTH2', 'onlyTH3']])
    objectNames = []
    allKeys = [k for k in inputDir.GetListOfKeys()]
    directoryKeys = [k for k in allKeys if r.TClass(k.GetClassName()).InheritsFrom(r.TDirectory.Class())]
    directoryNames = [k.GetName() for k in directoryKeys]
    def isTH1key(key) :
        kc = r.TClass(key.GetClassName())
        return kc.InheritsFrom(r.TH1.Class()) \
               and not kc.InheritsFrom(r.TH2.Class()) \
               and not kc.InheritsFrom(r.TH3.Class())
    def isTH2key(key) : return r.TClass(key.GetClassName()).InheritsFrom(r.TH2.Class())
    def isTH3key(key) : return r.TClass(key.GetClassName()).InheritsFrom(r.TH3.Class())
    if onlyTH1 : allKeys = [k for k in allKeys if isTH1key(k)]
    if onlyTH2 : allKeys = [k for k in allKeys if isTH2key(k)]
    if onlyTH3 : allKeys = [k for k in allKeys if isTH3key(k)]
    histoNames = [k.GetName() for k in allKeys if r.TClass(k.GetClassName()).InheritsFrom(r.TH1.Class())]
    if verbose : print histoNames
    for directory in directoryNames :
        histoNames += getAllHistoNames(inputDir.Get(directory), verbose, onlyTH1, onlyTH2, onlyTH3)
    return histoNames

def classifyHistoByName(histo, verbose=False) :
    cl = HistoNameClassifier(verbose)
    setattr(histo, 'type', cl.histoType(histo.GetName()))

def organizeHistosByType(histosByType = collections.defaultdict(list),
                         histosToOrganize = []) :
    "Fill a dictionary where the histos are keyed by HistoType"
    for h in histosToOrganize : histosByType[h.type].append(h)
    return histosByType
#
# testing
#
class KnownHistoTypes(unittest.TestCase) :
    def testMatchAllAvailabeAttrs(self):
        keys = ['pr', 'ch', 'var', 'syst']
        ht0 = HistoType(**dict(zip(keys, ('sr1','ee', 'l0_pt', 'NOM'))))
        ht1 = HistoType(**dict(zip(keys, ('sr1','ee', 'l0_pt', 'NOM'))))
        ht2 = HistoType(**dict(zip(keys, ('sr2','ee', 'l0_pt', 'NOM'))))
        ht3 = HistoType(**dict(zip(keys, (''   ,'ee', 'l0_pt', 'NOM'))))
        ht4 = HistoType(**dict(zip(keys, (''   ,'em', 'l0_pt', 'NOM'))))
        for (l,r), expRes in zip([(ht0,ht1), (ht0,ht2), (ht0,ht3), (ht3, ht1), (ht3, ht4)],
                                 [True,      False,     False,     True,       False]) :
            #print l.matchAllAvailabeAttrs(r)," <--> ",l," | ",r
            self.assertEqual(l.matchAllAvailabeAttrs(r), expRes)

class KnownHistoNames(unittest.TestCase) :
    def testAttrExtraction(self):
        "Verifiy that we are able to extract the parameters with weird histonames"
        hnc = HistoNameClassifier()
        for sr, c, v, sy in [('br4', 'em', 'dPhi_woSig_llmet_j','NOM'),
                             ('br4', 'em', 'dPhi_woSig_met_l0' ,'NOM'),
                             ('br4', 'em', 'dPhi_woSig_l0_l1'  ,'NOM'),
                             ('br4', 'em', 'l0_l1_pt'          ,'NOM'),
                             ('br4', 'em', 'll_M_fine'         ,'NOM'),
                             ('br4', 'em', 'll_M_fine.foo'     ,'NOM'),
                             ] :
            ht = hnc.histoType("%s_%s_%s_%s" % (sr, c, v, sy))
            self.assertEqual(all([l==r
                                  for l, r in zip([sr, c, v, sy],
                                                  [getattr(ht, a) for a in 'pr', 'ch', 'var', 'syst'])]),
                             True)

if __name__ == "__main__":
    unittest.main()

