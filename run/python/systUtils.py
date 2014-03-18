# Utility functions to handle systematic variations
#
# davide.gerbaudo@gmail.com
# March 2014

# todo : replace the various range(...) with rootUtils.getBinIndices()

import math
import numpy as np
from rootUtils import importRoot
r = importRoot()

from rootUtils import integralAndError

def fakeSystVariations() :
    "2x2x2=8 syst variations for the fake estimate, see DiLeptonMatrixMethod::systematic_names"
    return ['EL_RE_UP', 'EL_RE_DOWN', 'MU_RE_UP', 'MU_RE_DOWN',
            'EL_FR_UP', 'EL_FR_DOWN', 'MU_FR_UP', 'MU_FR_DOWN']
def mcObjectVariations() :
    "See definitions in SusyDefs.h:SusyNtSystNames, and active list in SusyPlotter::toggleStdSystematics()"
    return ['EES_Z_UP', 'EES_Z_DN',
            'EES_MAT_UP','EES_MAT_DN',
            'EES_PS_UP', 'EES_PS_DN',
            'EES_LOW_UP', 'EES_LOW_DN',
            'EER_UP', 'EER_DN',
            'MS_UP', 'MS_DN',
            'ID_UP', 'ID_DN',
            'JES_UP', 'JES_DN',
            'JER',
            'SCALEST_UP', 'SCALEST_DN',
            'RESOST',
            ]
def mcWeightVariations() :
    "See list at HftFiller::assignWeightVars()"
    return ['BKGMETHODUP' ,'BKGMETHODDOWN'
            ,'ETRIGREWUP' ,'ETRIGREWDOWN'
            ,'MTRIGREWUP' ,'MTRIGREWDOWN'
            ,'ESFUP'      ,'ESFDOWN'
            ,'MEFFUP'     ,'MEFFDOWN'
            ,'BJETUP'     ,'BJETDOWN'
            ,'CJETUP'     ,'CJETDOWN'
            ,'BMISTAGUP'  ,'BMISTAGDOWN'
            ,'XSUP'       ,'XSDOWN'
            ]
def mcWeightBranchname(mcWeightVariation='') : return 'syst_'+mcWeightVariation
def mcWeightBranches() : return [mcWeightBranchname(v) for v in mcWeightVariations()]

def getAllVariations() :
    return ['NOM'] + fakeSystVariations() + mcObjectVariations() + mcWeightVariations()

def fetchVariationHistos(input_fake_file=None, nominal_histo=None, variations=fakeSystVariations()) :
    nom_hname = nominal_histo.GetName()
    return dict([(v, input_fake_file.Get(nom_hname.replace('_NONE','_'+v))) for v in variations])
def computeFakeSysErr2(nominal_histo=None, vars_histos={}) :
    "Compute the bin-by-bin sum2 err including up&down fake systematic variations"
    def bc(h) : return [h.GetBinContent(b) for b in range(1, 1+h.GetNbinsX())]
    nom_bcs  = bc(nominal_histo)
    vars_bcs = dict([(v, bc(h)) for v, h in vars_histos.iteritems()])
    bins = range(nominal_histo.GetNbinsX())
    deltas = [[vars_bcs[v][b] - nom_bcs[b] for v in vars_bcs.keys()] for b in bins]
    def positive(ll) : return [l if not l<0.0 else 0.0 for l in ll ]
    def negative(ll) : return [l if     l<0.0 else 0.0 for l in ll ]
    def sumquad(ll) : return sum([l*l for l in ll])
    up_e2s = np.array([sumquad(positive(deltas[b])) for b in bins])
    do_e2s = np.array([sumquad(negative(deltas[b])) for b in bins])
    return {'up' : up_e2s, 'down' : do_e2s}
def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}
def computeFakeSysStatErr2(nominal_histo=None, vars_histos={}) :
    "Compute the bin-by-bin sum2 err including up&down fake systematic variations + stat. unc."
    print "refactor, use computeStatErr2 and computeSysErr2"
    def bc(h) : return [h.GetBinContent(b) for b in range(1, 1+h.GetNbinsX())]
    def be(h) : return [h.GetBinError(b) for b in range(1, 1+h.GetNbinsX())]
    nom_bcs  = bc(nominal_histo)
    nom_be2s = [e*e for e in be(nominal_histo)]
    vars_bcs = dict([(v, bc(h)) for v, h in vars_histos.iteritems()])
    bins = range(nominal_histo.GetNbinsX())
    deltas = [[vars_bcs[v][b] - nom_bcs[b] for v in vars_bcs.keys()] for b in bins]
    def positive(ll) : return [l if not l<0.0 else 0.0 for l in ll ]
    def negative(ll) : return [l if     l<0.0 else 0.0 for l in ll ]
    def sumquad(ll) : return sum([l*l for l in ll])
    up_e2s = np.array([sumquad(positive(deltas[b])) for b in bins]) + np.array(nom_be2s)
    do_e2s = np.array([sumquad(negative(deltas[b])) for b in bins]) + np.array(nom_be2s)
    return {'up' : up_e2s, 'down' : do_e2s}
def fetchFakeSysHistosAndComputeSysErr2(input_fake_file=None, nominal_histo=None) :
    vars_histos = fetchVariationHistos(input_fake_file, nominal_histo)
    return computeFakeSysStatErr2(nominal_histo, vars_histos)
def buildErrBandGraph(histo_tot_bkg, err2s) :
    h = histo_tot_bkg
    bins = range(1, 1+h.GetNbinsX())
    x = np.array([h.GetBinCenter (b) for b in bins])
    y = np.array([h.GetBinContent(b) for b in bins])
    ex_lo = ex_hi = np.array([0.5*h.GetBinWidth(b) for b in bins])
    ey_lo, ey_hi = np.sqrt(err2s['down']), np.sqrt(err2s['up'])
    gr = r.TGraphAsymmErrors(len(bins), x, y, ex_lo, ex_hi, ey_lo, ey_hi)
    gr.SetMarkerSize(0)
    gr.SetFillStyle(3004)
    gr.SetFillColor(r.kGray+3)
    gr.SetLineWidth(2)
    return gr
def buildErrBandRatioGraph(errband_graph) :
    gr = errband_graph.Clone()
    points = range(gr.GetN())
    xs     = np.array([gr.GetX()[i] for i in points])
    ys     = np.array([gr.GetY()[i] for i in points])
    eys_lo = np.array([abs(gr.GetErrorYlow (i)) for i in points])
    eys_hi = np.array([abs(gr.GetErrorYhigh(i)) for i in points])
    ys_lo  = ys - eys_lo
    ys_hi  = ys + eys_hi
    def absFracDeltaFromUnity(y_nom, y_var) : return abs(y_var/y_nom - 1.0) if y_nom else 0.0
    eys_lo = [absFracDeltaFromUnity(n, v) for n, v in zip(ys, ys_lo)]
    eys_hi = [absFracDeltaFromUnity(n, v) for n, v in zip(ys, ys_hi)]
    for p, x, ey_lo, ey_hi in zip(points, xs, eys_lo, eys_hi) :
        gr.SetPoint(p, x, 1.0) # TGraph does not have a SetPointY, so we need to set both x and y
        gr.SetPointEYlow (p, ey_lo)
        gr.SetPointEYhigh(p, ey_hi)
    return gr

#___________________________________________________________
# todo : move Group here (fake and simBkgs are Group objects)
def buildTotBackgroundHisto(histoFakeBkg=None, histosSimBkgs={}) :
    hTemplate = histoFakeBkg
    if not hTemplate :
        hTemplate = first(histosSimBkgs)
        print "warning, cannot use fake as template; histo name will be wrong, based on %s"%hTemplate.GetName()
    totBkg = hTemplate.Clone(hTemplate.GetName().replace('_fake','_totbkg'))
    totBkg.Reset()
    totBkg.Sumw2()
    allBackgrounds = dict((group, histo) for group, histo in [('fake',histoFakeBkg)]+[(g,h) for g,h in histosSimBkgs.iteritems()])
    for group, histo in allBackgrounds.iteritems() :
        if histo :
            totBkg.Add(histo)
            integral, error = integralAndError(histo)
            print "adding %.3f +/- %.3f  : %s (%s)" % (integral, error, group, histo.GetName())
        else : print "buildStatisticalErrorBand: missing %s"%group
    integral, error = integralAndError(totBkg)
    print "totBkg %.3f +/- %.3f  : %s" % (integral, error, totBkg.GetName())
    return totBkg
#___________________________________________________________
def buildStatisticalErrorBand(histoTotBkg= None) :
    return buildErrBandGraph(histoTotBkg, computeStatErr2(histoTotBkg))
#___________________________________________________________
def buildSystematicErrorBand(fake=None, nominalHistosSimBkg={}, variable='', selection='') :
    fake.setSystNominal()
    nominal_histo = buildTotBackgroundHisto(fake.getHistogram(variable, selection), nominalHistosSimBkg)
    vars_histos = dict((sys,
                        buildTotBackgroundHisto(fake.setSyst(sys).getHistogram(variable=variable, selection=selection),
                                                nominalHistosSimBkg))
                        for sys in fakeSystVariations())
    fakeErr2 = computeFakeSysErr2(nominal_histo=nominal_histo, vars_histos=vars_histos)
    return buildErrBandGraph(nominal_histo, fakeErr2)
def combineStatAndSystErrorBands(statErrBand, systErrBand) :
    sqrt = math.sqrt
    totErrBand = statErrBand.Clone() if statErrBand else systErrBand.Clone() if systErrBand else None
    if statErrBand and systErrBand :
        points = range(totErrBand.GetN())
        eys_stat_lo = np.array([abs(statErrBand.GetErrorYlow (i)) for i in points])
        eys_stat_hi = np.array([abs(statErrBand.GetErrorYhigh(i)) for i in points])
        eys_syst_lo = np.array([abs(systErrBand.GetErrorYlow (i)) for i in points])
        eys_syst_hi = np.array([abs(systErrBand.GetErrorYhigh(i)) for i in points])
        eys_lo = np.sqrt(eys_stat_lo*eys_stat_lo + eys_syst_lo*eys_syst_lo))
        eys_hi = np.sqrt(eys_stat_hi*eys_stat_hi + eys_syst_hi*eys_syst_hi))
        for p, ey_lo, ey_hi in zip(points, eys_lo, eys_hi) :
            totErrBand.SetPointEYlow (p, ey_lo)
            totErrBand.SetPointEYhigh(p, ey_hi)
    return totErrBand
