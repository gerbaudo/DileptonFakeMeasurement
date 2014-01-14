#!/bin/env python

# compute the variables needed to train a multivariate discriminant, and save them to a tree
#
# davide.gerbaudo@gmail.com
# Jan 2014

import array
import glob
import os
import ROOT as r
r.gROOT.SetStyle('Plain')
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
from math import cos, fabs, sin, sqrt, pi
from utils import first
from rootUtils import drawLegendWithDictKeys
from SampleUtils import colors
from kin import phi_mpi_pi, addTlv, computeMt, computeHt, computeMetRel, getDilepType

rootcoredir = os.environ['ROOTCOREDIR']
r.gROOT.LoadMacro(rootcoredir+'/scripts/load_packages.C+')
r.load_packages()

def buildFnamesDict(samples, dir='', tag='') :
    filenames = dict((s, glob.glob(dir+'/'+s+'_'+tag+'.root')) for s in samples)
    for s,v in filenames.iteritems() :
        if len(v)!=1 : print "skipping '%s', ambiguous filenames :\n%s"%(s, str(v))
    filenames = dict((s,v[0]) for s,v in filenames.iteritems() if len(v)==1 and v[0])
    return filenames

mt2 = r.mt2_bisect.mt2()
def computeMt2(a, b, met, zeroMass, lspMass) :
    pa    = array.array( 'd', [0.0 if zeroMass else a.M(), a.Px(), a.Py() ] )
    pb    = array.array( 'd', [0.0 if zeroMass else b.M(), b.Px(), b.Py() ] )
    pmiss = array.array( 'd', [0.0, met.Px(), met.Py() ] )
    mt2.set_momenta(pa, pb, pmiss)
    mt2.set_mn(lspMass)
    return mt2.get_mt2()
def computeMt2j(l0, l1, j0, j1, met, zeroMass=False, lspMass=0.0) :
    "As described in CMS-SUS-13-017"
    mt2_00 = computeMt2(l0+j0, l1+j1, met, zeroMass, lspMass)
    mt2_01 = computeMt2(l1+j0, l0+j1, met, zeroMass, lspMass)
    return min([mt2_00, mt2_01])
def computeMljj(l0, l1, j0, j1) :
    "Todo: extened combinatorics to N_j>2; good enough for now (we have few cases with >=3j)"
    jj = j0+j1
    dr0, dr1 = jj.DeltaR(l0), jj.DeltaR(l1)
    return (jj+l0).M() if dr0<dr1 else (jj+l1).M()

leafNames = ['pt0', 'pt1', 'mll', 'mtllmet', 'ht', 'metrel', 'dphill', 'detall', 'mt2j', 'mljj', 'dphijj', 'detajj']
structDecl  = 'struct vars { '
structDecl += ' '.join(["float %s;"%v for v in leafNames])
structDecl += " };"
r.gROOT.ProcessLine(structDecl)
vars = r.vars()
def resetVars(v) :
    for l in leafNames : setattr(v, l, 0.0)

def createOutTree(filenames, dilepChan, nJetChan, tag='', overwrite=False) :
    assert dilepChan in ['ee','mm','em']
    assert nJetChan in ['eq1j', 'ge2j']
    outFilenames = dict()
    for sample, filename in filenames.iteritems() :
        print sample
        outFilename = '/tmp/'+sample+'_'+dilepChan+'_'+nJetChan+'.root'
        if os.path.exists(outFilename) and not overwrite :
            outFilenames[sample] = outFilename
            continue
        outFile = r.TFile.Open(outFilename, 'recreate')
        outTree = r.TTree("training","Training tree")
        outTree.Branch('vars', vars, '/F:'.join(leafNames))
        outTree.SetDirectory(outFile)
        file = r.TFile.Open(filename)
        tree = file.Get(treename)
        print "processing %s (%d entries)"%(sample, tree.GetEntries())
        for iEvent, event in enumerate(tree) :
            resetVars(vars)
            l0 = addTlv(event.l0)
            l1 = addTlv(event.l1)
            met = addTlv(event.met)
            jets = [addTlv(j) for j in event.jets]
            lepts = [addTlv(l) for l in event.lepts]
            #lepts = filter(lambda l : lepIsSeparatedFromOther(l, [l0, l1]), lepts)
            pars = event.pars
            dilepType = getDilepType(event.l0, event.l1)
            nJets = len(jets)
            if dilepType != dilepChan : continue
            if nJets<1 or (nJets==1 and nJetChan is 'ge2j') : continue
            vars.pt0 = l0.p4.Pt()
            vars.pt1 = l1.p4.Pt()
            vars.mll = (l0.p4+l1.p4).M()
            vars.mtllmet = computeMt(l0.p4 + l1.p4, met.p4)
            vars.ht = computeHt(met.p4, [l0.p4, l1.p4]+[j.p4 for j in jets])
            vars.metrel = computeMetRel(met.p4, [l0.p4, l1.p4]+[j.p4 for j in jets])
            vars.dphill = fabs(phi_mpi_pi(l0.p4.DeltaPhi(l1.p4)))
            vars.detall = fabs(l0.p4.Eta() - l1.p4.Eta())
            # third lep, mljj,
            if nJets >1 :
                j0, j1 = jets[0], jets[1]
                vars.mt2j = computeMt2j(l0.p4, l1.p4, j0.p4, j1.p4, met.p4)
                vars.mljj = computeMljj(l0.p4, l1.p4, j0.p4, j1.p4)
                vars.dphijj = fabs(phi_mpi_pi(j0.p4.DeltaPhi(j1.p4)))
                vars.detajj = fabs(j0.p4.Eta() - j1.p4.Eta())
            outTree.Fill()
        print "filled ",outTree.GetEntries()," entries"
        outFile.Write()
        outFile.Close()
        outFilenames[sample] = outFile.GetName()
    return outFilenames


def train(sigFiles=[], bkgFiles=[], dilepChan='', nJetChan='') :
    assert dilepChan in ['ee','mm','em']
    assert nJetChan in ['eq1j', 'ge2j']
    bkgTree, sigTree = r.TChain('training'), r.TChain('training')
    for b in bkgFiles : bkgTree.Add(b)
    for s in sigFiles : sigTree.Add(s)
    if not bkgTree.GetEntries() or not sigTree.GetEntries() :
        print "TMVA cannot handle empty training trees...exiting"
        return
    r.TMVA.Tools.Instance()
    fileOut = r.TFile("testmva_%s_%s.root"%(dilepChan, nJetChan),"RECREATE")
    factory = r.TMVA.Factory("TMVAClassification", fileOut,
                             ":".join(["!V",
                                       "!Silent",
                                       "!Color",
                                       "!DrawProgressBar",
                                       "Transformations=I;D;P;G,D",
                                       "AnalysisType=Classification"]
                                      ))
    mvaVars = [l for l in leafNames if l not in ['pt0', 'pt1']]
    for v in mvaVars : factory.AddVariable(v, 'F')
    sigCut = r.TCut("")
    bkgCut = r.TCut("")
    factory.AddSignalTree    (sigTree)
    factory.AddBackgroundTree(bkgTree)
    factory.SetBackgroundWeightExpression("")
    factory.SetSignalWeightExpression("")
    factory.PrepareTrainingAndTestTree(sigCut, bkgCut,
                                       ":".join(["nTrain_Signal=0",
                                                 "nTrain_Background=0",
                                                 "SplitMode=Random",
                                                 "NormMode=NumEvents",
                                                 "V"
                                                 ]))
    method = factory.BookMethod(r.TMVA.Types.kBDT, "BDT",
                                ":".join(["!H",
                                          "V",
                                          "NTrees=200",
                                          "nEventsMin=50",
                                          "MaxDepth=3",
                                          "BoostType=AdaBoost",
                                          "AdaBoostBeta=0.5",
                                          "SeparationType=GiniIndex",
                                          "nCuts=20",
                                          "PruneMethod=NoPruning",
                                          ]))
    methodC = factory.BookMethod( r.TMVA.Types.kCuts, "Cuts",
                                  ":".join(["!H","!V","FitMethod=GA","EffSel"]))

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

# ------
# main
# todo: move to main func, add cmd-line opt, etc.
# ------

treename = 'SusySel'
tag = 'Jan_11' #'Jan_11'
basedir = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04_dev/SusyTest0/run/out/susysel/'

bkgSamples = ['wjets', 'zjets', 'ttbar', 'diboson', 'heavyflavor']
bkgSamples = ['ttbar']
sigSamples = ["Herwigpp_sM_wA_noslep_notauhad_WH_2Lep_%d"%d for d in range(1,1+1)]

bkgFilenanes = buildFnamesDict(bkgSamples, basedir+'/merged/', tag)
sigFilenanes = buildFnamesDict(sigSamples, basedir, tag)

for ll in ['ee','mm','em'] :
    for nj in ['eq1j', 'ge2j'] :
        print '-'*3+ll+nj+'-'*3
        bkgTrainFilenanes = createOutTree(bkgFilenanes, ll, nj, tag=tag)
        sigTrainFilenames = createOutTree(sigFilenanes, ll, nj, tag=tag)

        bkgFiles = bkgTrainFilenanes.values()
        sigFiles = [f for f in sigTrainFilenames.values()] # only one signal for now
        train(sigFiles, bkgFiles, ll, nj)






