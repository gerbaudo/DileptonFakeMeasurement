#!/bin/env python

# Utility classes and functions to manage samples and subsamples
#
# davide.gerbaudo@gmail.com
# Jan 2013

import glob, os, re, unittest
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

def guessSampleFromFilename(filename='', verbose=False) :
    if 'top_' in filename : return 'ttbar'
    elif 'Zjet_' in filename : return 'zjets'
    elif 'ZZ_' in filename \
         or 'WW_' in filename \
         or 'WZ_' in filename : return 'diboson'
    elif 'Wjet_' in filename : return 'wjets'
    elif 'physics_Egamma' in filename or 'physics_Muons' in filename : return 'data'
    elif 'WH_2Lep' in filename or 'WH_3Lep' in filename :
        match = re.search('(WH_\dLep_\d+)', filename)
        assert match, "failed to extract the WH pattern from %s"%filename
        return match.group()
    else :
        if verbose : print "cannot guess samplename for %s" % filename

def isDataSample(samplename) : return 'data' in samplename
def isSigSample(samplename) : return 'WH_' in samplename
def isBkgSample(samplename) : return not isDataSample(samplename) and not isSigSample(samplename)
def guessReqidFromFilename(filename='', verbose=False) :
    match = re.search('mc12\_8TeV\.(\d+)\.', filename)
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
        def valid(self) : return all([hasattr(self, a) for a in ModeAWhDbPar.fields]) and self.ds.isdigit()

    def __init__(self, filename=xsReaderDataDir()+'/modeA_WH_MC1eqMN2.txt') :
        self.entries = [e for e in [ModeAWhDbPar.Entry(l) for l in open(filename).readlines()] if e.valid()]
    def mc1Mn1ByReqid(self, reqid) :
        entry = next(e for e in self.entries if e.ds == reqid)
        return float(entry.mc1), float(entry.mn1)
    def allMc1(self) : return [e.mc1 for e in self.entries]
    def allMn1(self) : return [e.mn1 for e in self.entries]

class ModeAWhDbReqid :
    "Using the filelists, map reqids to samplenames"
    def __init__(self, filenames = []) :
        self.entries = {}
        if not filenames : filenames = glob.glob(basePathArea() + '/SusyTest0/run/filelist/wA_noslep_WH_*Lep*txt')
        for f in filenames :
            rootfile = open(f).read()
            reqid  = guessReqidFromFilename(rootfile)
            sample = guessSampleFromFilename(rootfile)
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
                       ,('176641', 'WH_3Lep_1')
                        ]
        db = ModeAWhDbReqid()
        for reqid, sample in knownValues :
            self.assertEqual(reqid, db.reqidBySample(sample))


if __name__ == "__main__":
    unittest.main()
