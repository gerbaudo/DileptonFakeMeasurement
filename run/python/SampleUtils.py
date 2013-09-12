#!/bin/env python

# Utility classes and functions to manage samples and subsamples
#
# davide.gerbaudo@gmail.com
# Jan 2013

import glob, os, re, unittest
from datasets import datasets, allGroups, activeDatasets
import ROOT as r
r.gROOT.SetBatch(1)

colors = {
    'ttbar'     : r.kRed+1,
    'zjets'     : r.kOrange-2,
    'wjets'     : r.kBlue-2,
    'diboson'   : r.kSpring+1,
    'singletop' : r.kAzure-4,
    'multijet'  : r.kWhite
    }

allGroups, datasets = allGroups(datasets), activeDatasets(datasets)

def guessGroupFromFilename(filename='') :
    "Guess group from filename, either merged or un-merged"
    group = next((g for g in allGroups if g+'_' in filename), None)
    # assume filelists are <...>+<dataset>.txt
    group = group if group else next((d.group
                                      for d in datasets
                                      if d.name+'.' in filename
                                      or d.name+'_' in filename), None)
    return group

def isDataSample(samplename) : return 'data' in samplename or 'period' in samplename
def isSigSample(samplename) : return 'WH_' in samplename
def isBkgSample(samplename) : return not isDataSample(samplename) and not isSigSample(samplename)
def guessReqidFromFilename(filename='', verbose=False) :
    match = re.search('mc12\_8TeV\.(\d+)\.', filename)
    if verbose and not match : print "'%s' does not contain mc12\_8TeV\.(\d+)\."%filename
    return match.group(1)


def basePathArea() :
    path = os.path.realpath(__file__)
    return path[:path.rfind('SusyTest0')]
def xsReaderDataDir(basePath='') :
    relPath = '/SusyXSReader/data'
    return (basePath if basePath else basePathArea()) + relPath

class ModeAWhDbPar :
    fields = ['ds', 'mc1', 'mn1', 'xsec', 'xsecSys']
    class Entry:
        def __init__(self, line) :
            """parse a line that is expected to be formatted as follow:
            DS MC1 MN1 xsec xsec_sys
            and store what is necessary.
            """
            line = line.strip()
            words = line.split()
            for a,w in zip(ModeAWhDbPar.fields, words) : setattr(self, a, w)
        def valid(self) :
            return all([hasattr(self, a) for a in ModeAWhDbPar.fields]) and self.ds.isdigit()

    def __init__(self, filename=xsReaderDataDir()+'/modeA_WH_MC1eqMN2.txt') :
        def isValidLine(l) :
            return len(l.strip()) and l.split()[0].isdigit()
        self.entries = [e for e in [ModeAWhDbPar.Entry(l)
                                    for l in open(filename).readlines() if isValidLine(l)]
                        if e.valid()]
    def mc1Mn1ByReqid(self, reqid) :
        entry = next(e for e in self.entries if e.ds == reqid)
        return float(entry.mc1), float(entry.mn1)
    def allMc1(self) : return [e.mc1 for e in self.entries]
    def allMn1(self) : return [e.mn1 for e in self.entries]

class ModeAWhDbReqid :
    "Using the filelists, map reqids to samplenames"
    def __init__(self, filenames = []) :
        self.entries = {}
        if not filenames :
            filenames = glob.glob(basePathArea()
                                  +'/SusyTest0/run/filelist/'
                                  +'Herwigpp_simplifiedModel_wA_noslep_WH_*Lep_*.txt')
        assert len(filenames),"cannot initialize ModeAWhDbReqid without filelists"
        for f in filenames :
            rootfile = open(f).read()
            if not rootfile :
                print "warning, emtpy filelist %s"%f
                continue
            reqid  = guessReqidFromFilename(rootfile)
            sample = guessGroupFromFilename(rootfile)
            if not sample : continue
            assert sample not in self.entries, "Cannot have several reqids with the same signal : %s, %s"%(sample, str([reqid, self.entries[sample]]))
            self.entries[sample] = reqid
    def reqidBySample(self, sample) :
        return self.entries[sample]
    def sampleByReqid(self, reqid) :
        return next((s for s,r in self.entries.iteritems() if r == reqid), None)

class ModeAWhDbMergedFake2Lreqid :
    "provide fake request ids to emulate merged samples in the 2L (mc1,mn1) plane"
    def __init__(self) :
        pass
    def reqidByMc1Mn1(self, mc1, mn1) :
        x, y = float(mc1), float(mn1)
        if   y<(-x+210) : return 1765700
        elif y<(-x+290) : return 1765800 if y>(x-210) else 1765900
        elif y<(-x+360) : return 1766000 if y>(x-210) else 1766100 if y>(x-290) else 1766200
        else            : return 1766300
        #if y<(x-100) else None # this is messing up also the bottom half? later on...

#
# testing
#
class KnownReqidModeAWhDb(unittest.TestCase) :
    def testMatchAllAvailabeAttrs(self) :
        knownValues = [ ('176574', (130.0 , 0.0))
                       ,('176575', (140.0, 10.0))
                        ]
        db = ModeAWhDbPar()
        for reqid, valuePair in knownValues :
            v1T, v2T = valuePair[0], valuePair[1]
            v1, v2  = db.mc1Mn1ByReqid(reqid)
            self.assertEqual(v1, v1T)
            self.assertEqual(v2, v2T)

class KnownEntriesModeAWhDbReqid(unittest.TestCase) :
    def testMatchAllAvailabeAttrs(self) :
        knownValues = [ ('176584', 'WH_2Lep_11')
                       ,('176581', 'WH_2Lep_8')
                        ]
        db = ModeAWhDbReqid()
        for reqid, sample in knownValues :
            self.assertEqual(reqid, db.reqidBySample(sample))


if __name__ == "__main__":
    unittest.main()
