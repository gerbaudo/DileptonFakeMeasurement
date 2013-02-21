#!/bin/env python

import collections
import ROOT as r
r.gROOT.SetBatch(1)
r.gInterpreter.GenerateDictionary('vector< vector< int > >', 'vector')

import gen
inputFname = '/tmp/gerbaudo/'\
             'mc12_8TeV.176575.Herwigpp_UEEE3_CTEQ6L1_simplifiedModel_wA_noslep_WH_2Lep_2.merge.NTUP_SUSY.e1702_s1581_s1586_r3658_r3549_p1328_tid01176849_00/'\
             'NTUP_SUSY.01176849._000001.root.1'
inputFname = '/tmp/gerbaudo/'\
             'mc12_8TeV.176584.Herwigpp_UEEE3_CTEQ6L1_simplifiedModel_wA_noslep_WH_2Lep_11.merge.NTUP_SUSY.e1702_s1581_s1586_r3658_r3549_p1328_tid01176858_00/'\
             'NTUP_SUSY.01176858._000001.root.1'

#inputFname = '/gdata/atlas/mrelich/HiggsinoTestSamples/Higgsino.000002.NTUP.root'

inFile = r.TFile.Open(inputFname)
tree = inFile.Get('susy')
#tree = inFile.Get('truth')

nEntries = tree.GetEntries()
print '%d entries' % nEntries
nEntriesToPrint = 20
printEveryFactor = 2
lastPrinted = 1

higgsPdg = 25
counter = collections.defaultdict(int)
printer = gen.particlePrinter()

for iEntry in xrange(nEntries) :
    tree.GetEntry(iEntry)
    pdg        = tree.mc_pdgId
    status     = tree.mc_status
    parents    = tree.mc_parent_index
    children   = tree.mc_child_index
#    for iH in range(len(pdg)) :
#        if pdg[iH] != higgsPdg : continue
#        print "higgs [%d] parents %s children %s"%(iH, [p for p in parents[iH]], [c for c in children[iH]])
        
        
    #iHiggs = next((i for i,p in enumerate(pdg) if p==higgsPdg), None)
    iHiggs = [i for i,p in enumerate(pdg) if p==higgsPdg]
    if not all([i<len(pdg) for i in iHiggs]) :
        print "invalid entry %d" %iEntry
        continue
    higgsChildren = [[pdg[i] if i<len(pdg) else 0 for i in hhc] for hhc in [children[i] for i in iHiggs] ]
    higgsParents = [[pdg[i] if i<len(pdg) else 0 for i in hhp] for hhp in [parents[i] for i in iHiggs] ]
    # there can be several intermediate higgs
    boringIhiggs = [ih for i, ih in enumerate(iHiggs) if len(higgsChildren[i])<2 or higgsPdg in higgsChildren[i]]
    interestingIhiggs = [i for i in iHiggs if i not in boringIhiggs]
    higgsChildren = [hc for hc,ih in zip(higgsChildren, iHiggs) if ih in interestingIhiggs]
    higgsParents = [hp for hp,ih in zip(higgsParents, iHiggs) if ih in interestingIhiggs]

    for hc in higgsChildren : counter[str(sorted(hc))] += 1
    if not interestingIhiggs : counter['nothing interesting'] += 1
    if len(interestingIhiggs) > 2 : counter['multiple interesting'] += 1
    if (iEntry/(lastPrinted if lastPrinted else 1) >= printEveryFactor) or iEntry<nEntriesToPrint:
        print '-'*4 + " event %d "%iEntry + '-'*4
        print "run: %d event %d"
        print "iHiggs: ",iHiggs
        print "interesting iHiggs: ",str(interestingIhiggs)
        print "higgsChildren: ",str(higgsChildren)
        print "higgsParents: ",str(higgsParents)
        lastPrinted = iEntry
    if iEntry < nEntriesToPrint : printer.uponAcceptance(tree)

cleanCounts = collections.defaultdict(int) #categorize them discarding the MC gen extra stuff
print "Counts for each decay"
for ch,count in sorted(counter.iteritems(), key=lambda (k,v): v, reverse=True):
    print "%s : %d" % (ch, count)
    children = eval(ch)
    if   any(p in children for p in [-24, +24])  : cleanCounts['WW']     += count
    elif any(p in children for p in [23])        : cleanCounts['ZZ']     += count
    elif any(p in children for p in [-15, +15])  : cleanCounts['tautau'] += count
    elif any(p in children for p in [-5, +5])    : cleanCounts['bbbar']  += count
    elif any(p in children for p in [-13, +13])  : cleanCounts['mumu']   += count
    else : print "unknown decay %s"%ch
totCategorized = sum(cleanCounts.values())
for ch,count in sorted(cleanCounts.iteritems(), key=lambda (k,v): v, reverse=True):
    print "%s : %d (%.2f%%)" % (ch, count, 100.*count/totCategorized)
