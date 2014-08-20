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

counter = collections.defaultdict(int)
printer = gen.particlePrinter()

for iEntry in xrange(nEntries) :
    tree.GetEntry(iEntry)
    pdg        = tree.mc_pdgId
    status     = tree.mc_status
    parents    = tree.mc_parent_index
    children   = tree.mc_child_index
    findHiggs = gen.findInterestingHiggsWithChiAndPar
    interestingIhiggs, higgsChildren, higgsParents = findHiggs(pdg, parents, children)
    for hc in higgsChildren : counter[str(sorted(hc))] += 1
    if not interestingIhiggs : counter['nothing interesting'] += 1
    if len(interestingIhiggs) > 2 : counter['multiple interesting'] += 1
    if (iEntry/(lastPrinted if lastPrinted else 1) >= printEveryFactor) or iEntry<nEntriesToPrint:
        print '-'*4 + " event %d "%iEntry + '-'*4
        print "run: %d event %d" % (tree.RunNumber, tree.EventNumber)
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
    label = gen.guessHdecayLabel(children)
    cleanCounts[label] += count
totCategorized = sum([v for k,v in cleanCounts.iteritems() if k!='unknown'])
for ch,count in sorted(cleanCounts.iteritems(), key=lambda (k,v): v, reverse=True):
    print "%s : %d (%.2f%%)" % (ch, count, 100.*count/totCategorized)
