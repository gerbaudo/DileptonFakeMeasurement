#!/bin/env python
import os, pickle, sys
import ROOT as r
r.gROOT.SetBatch(1)
r.gSystem.Load('/export/home/gerbaudo/workarea/Susy2013/SusyNtuple/StandAlone/libSusyNtuple.so')

ntMakerPklFname = 'ntMaker.pkl'
ntMakerLogFname = '/export/home/gerbaudo/workarea/SusyNt/subm_n0126/foo.log'
referenPklFname = 'referen.pkl'
referenLogFname = '/tmp/passBadMuonVeto.log'


class Event:
    def __init__(self, evtNum, badJet=False, goodVtx=False, badMuon=False) :
        for a in ['evtNum', 'badJet', 'goodVtx', 'badMuon'] : setattr(self, a, eval(a))
def parseNtMakerLog(filename) :
    print "parsing ",filename
    bjFlags = []
    vtFlags = []
    bmFlags = []
    for line in open(filename).readlines() :
        line = line.strip()
        badJet = 'BadJet evt :' in line
        goodVtx = 'GoodVtx evt :' in line
        badMuon = 'BadMuon evt :' in line
        if not badJet and not goodVtx and not badMuon : continue
        evtNum = int(line.split(':')[1])
        if badJet : bjFlags.append(Event(evtNum=evtNum, badJet = ('passed' in line)))
        elif goodVtx : vtFlags.append(Event(evtNum=evtNum, goodVtx = ('passed' in line)))
        elif badMuon : bmFlags.append(Event(evtNum=evtNum, badMuon = ('passed' in line)))
    print "number of flags : bj %d, vt %d, bm %d" % (len(bjFlags), len(vtFlags), len(bmFlags))
    commonEvents = set([j.evtNum for j in bjFlags])
    for ll in [vtFlags, bmFlags] : commonEvents.intersection_update([l.evtNum for l in ll])
    bjFlags = [j for j in bjFlags if j.evtNum in commonEvents]
    vtFlags = [v for v in vtFlags if v.evtNum in commonEvents]
    bmFlags = [m for m in bmFlags if m.evtNum in commonEvents]
    assert len(bjFlags)==len(vtFlags), "something went wrong with the parsing bj %d != vt %d" % (len(bjFlags), len(vtFlags))
    assert len(vtFlags)==len(bmFlags), "something went wrong with the parsing vt %d != bm %d" % (len(vtFlags), len(bmFlags))
    flags = [Event(j.evtNum, j.badJet, v.goodVtx, m.badMuon) for j,v,m in zip(bjFlags, vtFlags, bmFlags)]
    output = open(ntMakerPklFname, 'wb')
    pickle.dump(flags, output)
    output.close()
    return flags
def parseReferenceFile(filename) :
    print "parsing ",filename
    evtNumIndex = 2
    eventNumbers = [int(l.strip().split('*')[evtNumIndex])
                    for l in open(filename).readlines()]
    output = open(referenPklFname, 'wb')
    pickle.dump(eventNumbers, output)
    output.close()
    return eventNumbers

ntMakerFlags = pickle.load(open(ntMakerPklFname, 'rb')) if os.path.exists(ntMakerPklFname) else parseNtMakerLog(ntMakerLogFname)
refEvents =    pickle.load(open(referenPklFname, 'rb')) if os.path.exists(referenPklFname) else parseReferenceFile(referenLogFname)
print "%d events from ntMaker"%len(ntMakerFlags)
print "%d events from ref log"%len(refEvents)

assert len(ntMakerFlags)>len(refEvents),"expecting to have two extra events for ntMaker"
for e in ntMakerFlags:
    if e.evtNum not in refEvents :
        print "extra event:"
        print '\n'.join(["%s : %s"%(a, str(getattr(e,a))) for a in ['evtNum', 'badJet', 'goodVtx', 'badMuon']])


sys.exit()

#baseDir = '/gdata/atlas/ucintprod/SusyNt/mc12_n0127/user.sfarrell.mc12_8TeV.126893.Sherpa_CT10_lllnu_WZ.SusyNt.e1434_s1499_s1504_r3658_r3549_p1328_n0127/'
#filenames = ['user.sfarrell.102989._00001.susyNt.root'
#             ,'user.sfarrell.102989._00002.susyNt.root.1'
#             ,'user.sfarrell.102989._00003.susyNt.root'
#             ,'user.sfarrell.102989._00004.susyNt.root.1'
#             ,'user.sfarrell.102989._00005.susyNt.root.1'
#             ]
baseDir = '/export/home/gerbaudo/workarea/SusyNt/subm_n0126/'
filenames = ['mc12_8TeV.126893.Sherpa_CT10_lllnu_WZ_n0127.root',]

filenames = [baseDir+f for f in filenames]

nEntriesToPrint = 10
histo = None
for fn in filenames :
    f = r.TFile.Open(fn)
    h = f.Get('rawCutFlow')
    if histo : histo.Add(h)
    else : histo = h
    tree = f.Get('susyNt')
    for iEntry in xrange(nEntriesToPrint) :
        tree.GetEntry(iEntry)
        print "%d : %d %d" % (iEntry, tree.event.run, tree.event.event)
    f.ls()
    print "bin contents"
    print '\n    '.join(["%s %.1f (fail : %.1f) " % (h.GetXaxis().GetBinLabel(i),
                                              h.GetBinContent(i),
                                              h.GetBinContent(i)-h.GetBinContent(i-1) if i>1 else 0.)
                         for i in range(1, h.GetNbinsX()+1)])

print "total"
print '\n'.join(["%s %.1f (fail : %.1f) " % (histo.GetXaxis().GetBinLabel(i),
                                             histo.GetBinContent(i),
                                             histo.GetBinContent(i)-histo.GetBinContent(i-1) if i>1 else 0.)
                 for i in range(1, histo.GetNbinsX()+1)])


can = r.TCanvas()
can.cd()
histo.Draw()
can.Update()
can.SaveAs('cutflow.png')
