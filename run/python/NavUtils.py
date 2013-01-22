#!/bin/env/python

# Utility classes and function to navigate the histograms in a root file
#
# davide.gerbaudo@gmail.com
# Jan 2013

import collections, re
import ROOT as r

class HistoType(object):
    "Type of histogram, defined by plot region, channel, variable, syst"
    def __init__(self, pr='', ch='', var='', syst=''):
        for att in ['pr', 'ch', 'var', 'syst'] : setattr(self, att, eval(att))
    def sameas(self, rhs):
        return all([getattr(self,att)==getattr(rhs,att) for att in ['pr', 'ch', 'var', 'syst']])
    def __eq__(self, other) : return self.sameas(other)
    def __str__(self) : return ', '.join(["%s : %s"%(a, getattr(self,a)) for a in ['pr', 'ch', 'var', 'syst']])
    def __hash__(self) : return hash(self.__str__())


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
    "Extract PR+'_'+chan+'_'+name+'_'+sys from histoname and attach an HistoType attribute"
    n = histo.GetName()
    p = re.compile('(?P<pr>.*?)_'   # plot region (non greedy)
                   +'(?P<ch>.*?)_'  # channel (non greedy)
                   +'(?P<var>.*)_'  # var name (greedy, can contain '_')
                   +'(?P<syst>.*)') # last token, everything that's left
    match = p.search(n)
    if not match :
        if verbose : print "cannot classify %s" % n
        return
    kargs = dict([(g, match.group(g)) for g in ['pr', 'ch', 'var', 'syst']])
    setattr(histo, 'type', HistoType(**kargs))

def organizeHistosByType(histosByType = collections.defaultdict(list),
                         histosToOrganize = [], sampleName = '') :
    "Fill a dictionary where the histos are keyed by HistoType"
    for h in histosToOrganize :
        setattr(h, 'sample', sampleName)
        histosByType[h.type].append(h)
    return histosByType
