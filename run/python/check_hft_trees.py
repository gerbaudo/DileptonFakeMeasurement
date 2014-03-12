#!/bin/env python

# Sanity checks on the hft trees: yields and basic plots
#
# todo:
# - plot histos
# - blind sr and sr-within-pre
# - add data
# davide.gerbaudo@gmail.com
# March 2014

import datetime
import math
import optparse
import pprint
from rootUtils import (getMinMax
                       ,topRightLegend
                       ,importRoot
                       )
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)
from utils import first

from SampleUtils import colors
from CutflowTable import CutflowTable

def main():
    print "check_hft_trees.py"
    verbose = True
    groups = allGroups()
    selections = allRegions()
    variables = variablesToPlot()
    counters = bookCounters(groups, selections)
    histos = bookHistos(variables, groups, selections)
    samplesPerGroup = allSamplesAllGroups()
    for group, samplesGroup in samplesPerGroup.iteritems() :
        if verbose : print 1*' ',group
        hsGroup   = histos  [group]
        cntsGroup = counters[group]
        for sample in samplesGroup :
            if verbose : print 2*' ',sample
            file = r.TFile.Open(getFilename(group, sample))
            if not file or not file.IsOpen() :
                print "misssing '%s'"%getFilename(group, sample)
                continue
            tree = file.Get(getTreename(group, sample))
            if not tree :
                print "missing tree '%s' from '%s'"%(getTreename(group, sample), getFilename(group, sample))
                continue
            fillAndCount(hsGroup, cntsGroup, tree)
    if verbose : print 'done'
    countTotalBkg(counters)
    blindGroups  =  [g for g in counters.keys() if g!='data']
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

def dataFilename(samplename, inputdir='out/susyplot', syst='NOM') :
    return "%(d)s/%(sys)s_%(s)s.PhysCont.root"%{'d':inputdir, 's':samplename, 'sys':syst}
def dataTreename(samplename) :
    return "id_%(s)s.PhysCont" % {'s' : samplename}
def fakeFilename(samplename='periodX.physics_Astream', inputdir='out/fakepred/', syst='NOM') :
    return "%(d)s/%(sys)s_fake.%(s)s.PhysCont.root"%{'d':inputdir, 's':samplename, 'sys':syst}
def fakeTreename(samplename) :
    return "id_fake.%(s)s.PhysCont"%{'s':samplename}
def mcFilename(dsid=123456, inputdir='out/susyplot', syst='NOM') :
    return "%(d)s/%(sys)s_%(s)s.root"%{'d':inputdir, 's':dsid, 'sys':syst}
def mcTreename(dsid=123456) :
    return "id_%d"%dsid
def getFilename(group, sample) :
    if   group=='data' : return dataFilename(sample)
    elif group=='fake' : return fakeFilename(sample)
    else : return mcFilename(sample)
def getTreename(group, sample) :
    if   group=='data' : return dataTreename(sample)
    elif group=='fake' : return fakeTreename(sample)
    else : return mcTreename(sample)
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
def fillAndCount(histos, counters, tree) :
    selections = allRegions()
    selWeights = dict((s, r.TTreeFormula(s, selectionFormulas(s), tree)) for s in selections)
    for ev in tree :
        weight = tree.eventweight
        passSels = dict((s, selWeights[s].EvalInstance()) for s in selections)
        for s in selections :
            counters[s] += (weight if passSels[s] else 0.0)

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
    asg = dict([(k,v) for k,v in mcDatasetids().iteritems()]
               +[('data', dataSampleNames())]
               +[('fake', dataSampleNames())])
    asg = dict((k,v) for k,v in asg.iteritems() if k in allGroups())
    return asg
def variablesToPlot() :
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
    backgrounds = [g for g in allSamplesAllGroups().keys() if g!='signal' and g!='data']
    selections = first(counters).keys()
    counters['totBkg'] = dict((s, sum(counters[b][s] for b in backgrounds)) for s in selections)

if __name__=='__main__' :
    main()
