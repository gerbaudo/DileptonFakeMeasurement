#!/bin/env python

# Sanity checks on the hft trees: yields and basic plots
#
# todo:
# - fill syst histos
# - combine errors
# davide.gerbaudo@gmail.com
# March 2014

import datetime
import math
import optparse
import os
import pprint
from rootUtils import (getMinMax
                       ,topRightLegend
                       ,importRoot
                       ,setWhPlotStyle
                       )
r = importRoot()
from utils import (first
                   ,mkdirIfNeeded
                   ,filterWithRegexp
                   )

import SampleUtils
from CutflowTable import CutflowTable
import systUtils

usage="""

This code is used either (1) to fill the histos, or (2) to make plots
and tables. The output of (1) is used as input of (2).

Required inputs: HFT trees produced with `submitJobs.py ... --with-hft`
See cmd/hft.txt for details.

Example usage ('fill' mode):
%prog \\
 --syst EES_Z_UP
 --input-gen  out/susyplot/merged/ \\
 --input-fake out/fakepred/merged/ \\
 --output_dir     out/hft/             \\
 --verbose
 >& log/hft/check_hft_trees_fill.log

Example usage ('plot' mode):
%prog \\
 --syst ANY            \\
 --input-dir out/hft/  \\
 --verbose             \\
 >& log/hft/check_hft_trees_plot.log"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-f', '--input-fake', help='location hft trees for fake')
    parser.add_option('-g', '--input-gen', help='location hft trees for everything else')
    parser.add_option('-i', '--input-dir')
    parser.add_option('-o', '--output-dir')
    parser.add_option('-s', '--syst', help="variations to process (default all). Give a list or say 'weight', 'object', or 'fake'")
    parser.add_option('-e', '--exclude', help="skip some systematics, example 'EL_FR_.*'")
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    parser.add_option('-l', '--list-systematics', action='store_true', default=False, help='list what is already in output_dir')
    parser.add_option('-L', '--list-all-systematics', action='store_true', default=False, help='list all possible systematics')

    (opts, args) = parser.parse_args()
    if opts.list_all_systematics :
        print "All systematics:\n\t%s"%'\n\t'.join(systUtils.getAllVariations())
        return
    if opts.list_systematics :
        print listExistingSyst(opts.input_dir)
        return
    inGenSpecified, inDirSpecified = opts.input_gen!=None, opts.input_dir!=None
    eitherMode = inGenSpecified != inDirSpecified
    if not eitherMode : parser.error("Run either in 'fill' or 'plot' mode")
    mode = 'fill' if inGenSpecified else 'plot' if inDirSpecified else None
    requiredOptions = (['input_fake', 'input_gen', 'output_dir'] if mode=='fill' else ['input_dir', 'output_dir'])
    allOptions = [x.dest for x in parser._get_all_options()[1:]]
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) :
        parser.error('Missing required option\n'+'\n'.join(["%s : %s"%(o, getattr(opts, o)) for o in requiredOptions]))
    if opts.verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    if   mode=='fill' : runFill(opts)
    elif mode=='plot' : runPlot(opts)

def runFill(opts) :
    inputFakeDir = opts.input_fake
    inputGenDir  = opts.input_gen
    outputDir    = opts.output_dir
    sysOption    = opts.syst
    excludedSyst = opts.exclude
    verbose      = opts.verbose

    mkdirIfNeeded(outputDir)
    systematics = ['NOM']
    anySys = sysOption==None
    if   sysOption=='fake'   or anySys : systematics += systUtils.fakeSystVariations()
    elif sysOption=='object' or anySys : systematics += systUtils.mcObjectVariations()
    elif sysOption=='weight' or anySys : systematics += systUtils.mcWeightVariations()
    elif sysOption in systUtils.getAllVariations() : systematics = [sysOption]
    else : raise ValueError("Invalid syst %s"%sysOption)
    if excludedSyst : systematics = [s for s in systematics if s not in filterWithRegexp(systematics, excludedSyst)]

    for syst in systematics :
        if verbose : print '---- filling ',syst
        samplesPerGroup = allSamplesAllGroups()
        [s.setSyst(syst) for g, samples in samplesPerGroup.iteritems() for s in samples]
        counters, histos = countAndFillHistos(samplesPerGroup=samplesPerGroup, syst=syst, verbose=verbose, outdir=outputDir)
        printCounters(counters)
        saveHistos(samplesPerGroup, histos, outputDir, verbose)

def runPlot(opts) :
    inputDir   = opts.input_dir
    outputDir  = opts.output_dir

    mkdirIfNeeded(outputDir)
    histos = fetchHistos(inputDir)
    plotHistos(histos, outputDir)

def countAndFillHistos(samplesPerGroup={}, syst='', verbose=False, outdir='./') :

    selections = allRegions()
    variables = variablesToPlot()

    mcGroups, fakeGroups = mcDatasetids().keys(), ['fake']
    objVariations, weightVariations, fakeVariations = systUtils.mcObjectVariations(), systUtils.mcWeightVariations(), systUtils.fakeSystVariations()
    def groupIsRelevantForSys(g, s) :
        isRelevant = (s=='NOM' or (g in mcGroups and s in objVariations+weightVariations) or (g in fakeGroups and s in fakeVariations))
        if verbose and not isRelevant : print "skipping %s for %s"%(g, s)
        return isRelevant
    def dropIrrelevantGroupsForThisSys(groups, sys) : return dict((g, samples) for g, samples in groups.iteritems() if groupIsRelevantForSys(g, syst))
    samplesPerGroup = dropIrrelevantGroupsForThisSys(samplesPerGroup, syst)
    def dropSamplesWithoutTree(samples) : return [s for s in samples if s.hasInputHftTree(msg='Warning! ')]
    samplesPerGroup = dict((g, dropSamplesWithoutTree(samples)) for g, samples in samplesPerGroup.iteritems())
    def dropGroupsWithoutSamples(groups) : return dict((g, samples) for g, samples in groups.iteritems() if len(samples))
    samplesPerGroup = dropGroupsWithoutSamples(samplesPerGroup)

    groups = samplesPerGroup.keys()
    counters = bookCounters(groups, selections)
    histos = bookHistos(variables, groups, selections)
    for group, samplesGroup in samplesPerGroup.iteritems() :
        if verbose : print 1*' ',group
        hsGroup   = histos  [group]
        cntsGroup = counters[group]
        for sample in samplesGroup :
            if verbose : print 2*' ',sample.name
            fillAndCount(hsGroup, cntsGroup, sample)
    if verbose : print 'done'
    return counters, histos

def printCounters(counters):
    countTotalBkg(counters)
    blindGroups   = [g for g in counters.keys() if g!='data']
    unblindGroups = [g for g in counters.keys()]
    tableSr  = CutflowTable(samples=blindGroups,   selections=signalRegions(), countsSampleSel=counters)
    tablePre = CutflowTable(samples=blindGroups,   selections=controlRegions(), countsSampleSel=counters)
    tableBld = CutflowTable(samples=unblindGroups, selections=blindRegions(), countsSampleSel=counters)
    for table in [tableSr, tablePre, tableBld] : table.nDecimal = 6
    print 4*'-',' sig regions ',4*'-'
    print tableSr.csv()
    print 4*'-',' pre regions ',4*'-'
    print tablePre.csv()
    print 4*'-',' blind regions ',4*'-'
    print tableBld.csv()
def plotHistos(histos, plotdir='./') :
    blindGroups = [g for g in histos.keys() if g!='data']
    countersFromHist = dict((g, dict((s, histos[g][s]['mljj'].Integral()) for s in controlRegions())) for g in blindGroups)
    tableHist = CutflowTable(samples=blindGroups, selections=controlRegions(), countsSampleSel=countersFromHist)
    tableHist.nDecimal = 6
    print 4*'-',' mljj histograms ',4*'-'
    print tableHist.csv()
    histosPerSel = dict((s, dict((g, histos[g][s]['mljj']) for g in histos.keys())) for s in first(histos).keys())
    for s, hs in histosPerSel.iteritems() : plotHistos(s, hs, plotdir)

#___________________________________________________________
class BaseSampleGroup(object) :
    def __init__(self, name) :
        self.name = name
        self.setSyst()
    @property
    def label(self) : return self.groupname if hasattr(self, 'groupname') else self.name
    @property
    def isFake(self) : return self.label=='fake'
    @property
    def isData(self) : return self.label=='data'
    @property
    def isMc(self) : return not (self.isFake or self.isData)
    def isNeededForSys(self, sys) :
        return (sys=='NOM'
                or (self.isMc and sys in systUtils.mcWeightVariations())
                or (self.isMc and sys in systUtils.mcObjectVariations())
                or (self.isFake and sys in systUtils.fakeSystVariations()))
    def setSyst(self, sys='NOM') :
        nominal = 'NOM' # do we have differnt names for nom (mc vs fake)?
        self.isObjSys    = sys in systUtils.mcObjectVariations()
        self.isWeightSys = sys in systUtils.mcWeightVariations()
        self.isFakeSys   = sys in systUtils.fakeSystVariations()
        def nameObjectSys(s) : return s if self.isMc else nominal
        def nameWeightSys(s) : return s if self.isMc else nominal
        def nameFakeSys(s) : return s if self.isFake else nominal
        def identity(s) : return s
        sysNameFunc = nameObjectSys if self.isObjSys else nameWeightSys if self.isWeightSys else nameFakeSys if self.isFakeSys else identity
        self.syst = sysNameFunc(sys)
        return self
#___________________________________________________________
class Sample(BaseSampleGroup) :
    def __init__(self, name, groupname) :
        super(Sample, self).__init__(name) # this is either the name (for data and fake) or the dsid (for mc)
        self.groupname = groupname
        self.setHftInputDir()
    def setHftInputDir(self, dir='') :
        useDefaults = not dir
        defaultDir = 'out/fakepred' if self.isFake else 'out/susyplot'
        self.hftInputDir = defaultDir if useDefaults else dir
        return self
    @property
    def weightLeafname(self) : return 'eventweight' if not self.isWeightSys else systUtils.mcWeightBranchname(self.syst)
    @property
    def filenameHftTree(self) :
        def dataFilename(sample, dir, sys) : return "%(dir)s/%(sys)s_%(sam)s.PhysCont.root" % {'dir':dir, 'sam':sample, 'sys':sys}
        def fakeFilename(sample, dir, sys) : return "%(dir)s/%(sys)s_fake.%(sam)s.PhysCont.root" % {'dir':dir, 'sam':sample, 'sys':sys}
        def mcFilename  (sample, dir, sys) : return "%(dir)s/%(sys)s_%(dsid)s.root" % {'dir':dir, 'sys':sys, 'dsid':sample}
        fnameFunc = dataFilename if self.isData else fakeFilename if self.isFake else mcFilename
        sys = self.syst if self.isMc and self.isObjSys else 'NOM'
        return fnameFunc(self.name, self.hftInputDir, sys)
    @property
    def hftTreename(self) :
        def dataTreename(samplename) : return "id_%(s)s.PhysCont" % {'s' : samplename}
        def fakeTreename(samplename) : return "id_fake.%(s)s.PhysCont"%{'s':samplename}
        def mcTreename(dsid=123456) :  return "id_%d"%dsid
        getTreename = dataTreename if self.isData else fakeTreename if self.isFake else mcTreename
        return getTreename(self.name)
    def hasInputHftFile(self, msg) :
        filename = self.filenameHftTree
        isThere = os.path.exists(filename)
        if not isThere : print msg+"%s %s missing : %s"%(self.groupname, self.name, filename)
        return isThere
    def hasInputHftTree(self, msg='') :
        treeIsThere = False
        if self.hasInputHftFile(msg) :
            filename, treename = self.filenameHftTree, self.hftTreename
            inputFile = r.TFile.Open(filename) if self.hasInputHftFile(msg) else None
            if inputFile :
                if inputFile.Get(treename) : treeIsThere = True
                else : print msg+"%s %s missing tree '%s' from %s"%(self.groupname, self.name, treename, filename)
            inputFile.Close()
        return treeIsThere
    def group(self) :
        return Group(self.groupname).setSyst(self.syst)
#___________________________________________________________
class Group(BaseSampleGroup) :
    def __init__(self, name) :
        super(Group, self).__init__(name)
        self.setSyst()
        self.setHistosDir()
    def setHistosDir(self, dir='') :
        self.histosDir = dir if dir else 'out/hft'
        return self
    @property
    def filenameHisto(self) :
        def dataFilename(group, dir, sys) : return "%(dir)s/%(sys)s_%(gr)s.PhysCont.root" % {'dir':dir, 'gr':group, 'sys':sys}
        def fakeFilename(group, dir, sys) : return "%(dir)s/%(sys)s_fake.%(gr)s.PhysCont.root" % {'dir':dir, 'gr':group, 'sys':sys}
        def mcFilename  (group, dir, sys) : return "%(dir)s/%(sys)s_%(gr)s.root" % {'dir':dir, 'sys':sys, 'gr':group}
        fnameFunc = dataFilename if self.isData else fakeFilename if self.isFake else mcFilename
        return fnameFunc(self.name, self.histosDir, self.syst)

#___________________________________________________________
def allGroups(noData=False, noSignal=True) :
    return ([k for k in mcDatasetids().keys() if k!='signal' or not noSignal]
            + ([] if noData else ['data'])
            + ['fake']
            )

def llPairs() : return ['ee', 'em', 'mm']
def njetSelections() : return ['1jet', '23jets']
def signalRegions() :
    return ["%(ll)sSR%(nj)s"%{'ll':ll, 'nj':nj} for ll in llPairs() for nj in njetSelections()]
def controlRegions() :
    return ['pre'+r for r in signalRegions()]
def blindRegions() :
    return ['bld'+r for r in signalRegions()]
def allRegions() :
    return signalRegions() + controlRegions() + blindRegions()
def selectionFormulas(sel) :
    ee, em, mm = 'isEE', 'isEMU', 'isMUMU'
    pt32  = '(lept1Pt>30000.0 && lept2Pt>20000.0)'
    pt33  = '(lept1Pt>30000.0 && lept2Pt>30000.0)'
    j1    = '(L2nCentralLightJets==1)'
    j23   = '(L2nCentralLightJets==2 || L2nCentralLightJets==3)'
    vetoZ = '(L2Mll<(91200.0-10000.0) || L2Mll>(91200.0+10000.))'
    dEll  = '(TMath::Abs(deltaEtaLl)<1.5)'
    ss    = '(!isOS || L2qFlipWeight!=1.0)' # ssOrQflip
    mlj1  = 'mlj < 90000.0'
    mlj2  = 'mljj<120000.0'
    formulas = {
        'eeSR1jet'   : '('+ee+' && '+ss+' && '+j1 +' && '+pt32+' && '+vetoZ+' && '+mlj1+' && L2METrel>55000.0 && Ht>200000.0)',
        'eeSR23jets' : '('+ee+' && '+ss+' && '+j23+' && '+pt32+' && '+vetoZ+' && '+mlj2+' && L2METrel>30000.0 &&                mtmax>100000.0)',
        'mmSR1jet'   : '('+mm+' && '+ss+' && '+j1 +' && '+pt32+' && '+dEll +' && '+mlj1+' &&                     Ht>200000.0 && mtmax>100000.0)',
        'mmSR23jets' : '('+mm+' && '+ss+' && '+j23+' && '+pt33+' && '+dEll +' && '+mlj2+' &&                     Ht>220000.0)',
        'emSR1jet'   : '('+em+' && '+ss+' && '+j1 +' && '+pt33+' && '+dEll +' && '+mlj1+' && mtllmet>120000.0 &&                mtmax>110000.0)',
        'emSR23jets' : '('+em+' && '+ss+' && '+j23+' && '+pt33+' && '+dEll +' && '+mlj2+' && mtllmet>110000.0 )',
        }
    for f in formulas.keys() :
        formulas['pre'+f] = formulas[f].replace(mlj1, '1').replace(mlj2, '1')
    mlj1Not, mlj2Not = mlj1.replace('<','>'), mlj2.replace('<','>')
    for f in formulas.keys() :
        formulas['bld'+f] = formulas[f].replace(mlj1, mlj1Not).replace(mlj2, mlj2Not)
    return formulas[sel]

def fillAndCount(histos, counters, sample, blind=True) :
    group    = sample.group
    filename = sample.filenameHftTree
    treename = sample.hftTreename
    file = r.TFile.Open(filename)
    tree = file.Get(treename)
    selections = allRegions()
    selWeights = dict((s, r.TTreeFormula(s, selectionFormulas(s), tree)) for s in selections)
    weightFormula = r.TTreeFormula('weightFormula', sample.weightLeafname, tree)
    for ev in tree :
        weight = weightFormula.EvalInstance()
        passSels = dict((s, selWeights[s].EvalInstance()) for s in selections)
        for s in selections : counters[s] += (weight if passSels[s] else 0.0)
        for sr, pr in zip(signalRegions(), controlRegions()) :
            passPre, passSig = passSels[pr], passSels[sr]
            fillHisto = passPre and not passSig if (group=='data' and blind) else passPre
            oneJet = tree.L2nCentralLightJets==1
            mev2gev = 1.0e-3
            mljj = mev2gev*(tree.mlj if oneJet else tree.mljj)
            if fillHisto : histos[pr]['mljj'].Fill(mljj, weight)
    file.Close()

def mcSystematics() :
    return ['NOM', 'EER_DN', 'EER_UP', 'EES_LOW_DN', 'EES_LOW_UP',
            'EES_MAT_DN', 'EES_MAT_UP', 'EES_PS_DN', 'EES_PS_UP',
            'EES_Z_DN', 'EES_Z_UP', 'ID_DN', 'ID_UP', 'JER', 'JES_DN',
            'JES_UP', 'MS_DN', 'MS_UP', 'RESOST', 'SCALEST_DN',
            'SCALEST_UP',]
def fakeSystematics() :
    return ['NOM', 'EL_FR_DOWN', 'EL_FR_UP', 'EL_RE_DOWN', 'EL_RE_UP',
            'MU_FR_DOWN', 'MU_FR_UP', 'MU_RE_DOWN', 'MU_RE_UP',]
def dataSampleNames() :
    return ["period%(period)s.physics_%(stream)s"%{'period':p, 'stream':s}
            for p in ['A','B','C','D','E','G','H','I','J','L']
            for s in ['Egamma','Muons']]
def mcDatasetids() :
    "encode the grouping we use to make HFT plots; from HistFitter_TreeCreator.py"
    return {
        'Higgs' : [160155, 160205, 160255, 160305, 160505, 160555,
                   160655, 160705, 160755, 160805, 161005, 161055,
                   161105, 161155, 161305, 161555, 161566, 161577,
                   161675, 161686, 161697, 161708, 161719, 161730,
                   161805, 167418, 169072],
        'Top' : [110001, 108346, 119353, 174830, 174831, 119355,
                   174832, 174833, 179991, 179992, 119583, 169704,
                   169705, 169706, 158344],
        'WW' : [177997, 183734, 183736, 183738, 169471, 169472,
                169473, 169474, 169475, 169476, 169477, 169478,
                169479, 126988, 126989, 167006],
        'ZV' : [179974, 179975, 183585, 183587, 183589, 183591,
                183735, 183737, 183739, 167007, 177999, 183586,
                183588, 183590, 126894, 179396, 167008],
        'Zjets' : [178354, 178355, 178356, 178357, 178358, 178359,
                   178360, 178361, 178362, 178363, 178364, 178365,
                   178366, 178367, 178368, 178369, 178370, 178371,
                   178372, 178373, 178374, 178375, 178376, 178377,
                   178378, 178379, 178380, 178381, 178382, 178383,
                   178384, 178385, 178386, 178387, 178388, 178389,
                   178390, 178391, 178392, 178393, 178394, 178395,
                   178396, 178397, 178398, 178399, 178400, 178401,
                   178402, 178403, 178404, 178405, 178406, 178407,
                   117650, 117651, 117652, 117653, 117654, 117655,
                   110805, 110806, 110807, 110808, 110817, 110818,
                   110819, 110820, 117660, 117661, 117662, 117663,
                   117664, 117665, 110809, 110810, 110811, 110812,
                   110821, 110822, 110823, 110824, 117670, 117671,
                   117672, 117673, 117674, 117675, 110813, 110814,
                   110815, 110816, 110825, 110826, 110827, 110828],
        'signal' : [177501, 177502, 177503, 177504, 177505, 177506,
                    177507, 177508, 177509, 177510, 177511, 177512,
                    177513, 177514, 177515, 177516, 177517, 177518,
                    177519, 177520, 177521, 177522, 177523, 177524,
                    177525, 177526]
        }
def allSamplesAllGroups() :
    asg = dict( [(k, [Sample(groupname=k, name=d) for d in v]) for k,v in mcDatasetids().iteritems()]
               +[('data', [Sample(groupname='data', name=s) for s in dataSampleNames()])]
               +[('fake', [Sample(groupname='fake', name=s) for s in dataSampleNames()])])
    asg = dict((k,v) for k,v in asg.iteritems() if k in allGroups())
    return asg
def stackedGroups() :
    return [g for g in allSamplesAllGroups().keys() if g not in ['data', 'signal']]

def variablesToPlot() :
    return ['mljj']
    return ['pt0','pt1','mll','mtmin','mtmax','mtllmet','ht','metrel','dphill','detall',
            'mt2j','mljj','dphijj','detajj']
def histoSuffix(sample, selection) : return "%s_%s"%(sample, selection)
def bookHistos(variables, samples, selections) :
    "book a dict of histograms with keys [sample][selection][var]"
    def histo(variable, suffix) :
        s = suffix
        twopi = +2.0*math.pi
        h = None
        if   v=='pt0'     : h = r.TH1F('h_pt0_'    +s, ';p_{T,l0} [GeV]; entries/bin',          25, 0.0, 250.0)
        elif v=='pt1'     : h = r.TH1F('h_pt1_'    +s, ';p_{T,l1} [GeV]; entries/bin',          25, 0.0, 250.0)
        elif v=='mll'     : h = r.TH1F('h_mll_'    +s, ';m_{l0,l1} [GeV]; entries/bin',         25, 0.0, 250.0)
        elif v=='mtmin'   : h = r.TH1F('h_mtmin_'  +s, ';m_{T,min}(l, MET) [GeV]; entries/bin', 25, 0.0, 400.0)
        elif v=='mtmax'   : h = r.TH1F('h_mtmax_'  +s, ';m_{T,max}(l, MET) [GeV]; entries/bin', 25, 0.0, 400.0)
        elif v=='mtllmet' : h = r.TH1F('h_mtllmet_'+s, ';m_{T}(l+l, MET) [GeV]; entries/bin',   25, 0.0, 600.0)
        elif v=='ht'      : h = r.TH1F('h_ht_'     +s, ';H_{T} [GeV]; entries/bin',             25, 0.0, 800.0)
        elif v=='metrel'  : h = r.TH1F('h_metrel_' +s, ';MET_{rel} [GeV]; entries/bin',         25, 0.0, 300.0)
        elif v=='dphill'  : h = r.TH1F('h_dphill_' +s, ';#Delta#phi(l, l) [rad]; entries/bin',  25, 0.0, twopi)
        elif v=='detall'  : h = r.TH1F('h_detall_' +s, ';#Delta#eta(l, l); entries/bin',        25, 0.0, +3.0 )
        elif v=='mt2j'    : h = r.TH1F('h_mt2j_'   +s, ';m^{J}_{T2} [GeV]; entries/bin',        25, 0.0, 500.0)
        elif v=='mljj'    : h = r.TH1F('h_mljj_'   +s, ';m_{ljj} [GeV]; entries/bin',           25, 0.0, 500.0)
        elif v=='dphijj'  : h = r.TH1F('h_dphijj_' +s, ';#Delta#phi(j, j) [rad]; entries/bin',  25, 0.0, twopi)
        elif v=='detajj'  : h = r.TH1F('h_detajj_' +s, ';#Delta#eta(j, j); entries/bin',        25, 0.0, +3.0 )
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        return h
    return dict([(s,
                  dict([(sel,
                         dict([(v, histo(v, histoSuffix(s, sel))) for v in variables]))
                         for sel in selections]))
                 for s in samples])
def bookCounters(samples, selections) :
    "book a dict of counters with keys [sample][selection]"
    return dict((s, dict((sel, 0.0) for sel in selections)) for s in samples)
def countTotalBkg(counters={'sample' : {'sel':0.0}}) :
    backgrounds = [g for g in counters.keys() if g!='signal' and g!='data']
    selections = first(counters).keys()
    counters['totBkg'] = dict((s, sum(counters[b][s] for b in backgrounds)) for s in selections)
def getGroupColor(g) :
    oldColors = SampleUtils.colors
    colors = dict((g,c) for g,c in [(k,v) for k,v in oldColors.iteritems()] + [('signal',r.kRed), ('WW',r.kBlue), ('Higgs',r.kYellow)])
    def hftGroup2stdGroup(_) :
        fromTo = {'Zjets':'zjets', 'Top':'ttbar', 'ZV':'diboson',}
        return fromTo[_] if _ in fromTo else _
    g = hftGroup2stdGroup(g)
    return colors[g]

def plotHistos(selection='', histos={}, outdir='./') :
    setWhPlotStyle()
    padMaster = first(histos)
    can = r.TCanvas('can_'+padMaster.GetName(), padMaster.GetTitle(), 800, 600)
    can.cd()
    can._hists = [padMaster]
    padMaster.Draw('axis')
    stack = r.THStack('stack_'+padMaster.GetName(), '')
    can._hists.append(stack)
    leg = topRightLegend(can, 0.275, 0.475, shift=-0.025)
    can._leg = leg
    leg.SetBorderSize(0)
#     leg.AddEntry(h_data, dataSample(), 'P')
#     leg.AddEntry(h_bkg, 'sm', 'L')
    for g in stackedGroups() :
        if g not in histos : continue
        h = histos[g]
        h.SetFillColor(getGroupColor(g))
        h.SetLineColor(h.GetFillColor())
        stack.Add(h)
        can._hists.append(h)
    for g in stackedGroups()[::-1] : leg.AddEntry(histos[g], g, 'F') # stack goes b-t, legend goes t-b
    stack.Draw('hist')
#     leg.AddEntry(err_band, 'Uncertainty', 'F')
    leg.Draw('same')
    can.Update() # force stack to create padMaster
    hStack = stack.GetHistogram()
    padMaster.SetMaximum(1.1*max([h.GetMaximum() for h in [hStack]+histos.values()]))
    can.Update()
    for ext in ['png'] : can.SaveAs(outdir+'/'+can.GetName()+'.'+ext)

def listExistingSyst(dir) :
    print "listing systematics from ",dir
    print "...not implemented..."
def saveHistos(samplesPerGroup={}, histosPerGroup={}, outdir='./', verbose=False) :
    for groupname, histos in histosPerGroup.iteritems() :
        group = first(samplesPerGroup[groupname]).group().setHistosDir(outdir)
        outFilename = group.filenameHisto
        if verbose : print "creating file %s"%outFilename
        file = r.TFile.Open(outFilename, 'recreate')
        file.cd()
        for g, histosPerSel in histosPerGroup.iteritems() :
            for sel, histos in histosPerSel.iteritems() :
                for var, h in histos.iteritems() : h.Write()
        file.Close()

def fetchHistos(dir) :
    return "to be implemented"

if __name__=='__main__' :
    main()
