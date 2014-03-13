# Utility functions to handle systematic variations
#
# davide.gerbaudo@gmail.com
# March 2014

# todo : replace the various range(...) with rootUtils.getBinIndices()

import numpy as np
from rootUtils import importRoot
r = importRoot()

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
    "See list at SusyPlotter::computeWeightVariations()"
    return ['qflipUp', 'qflipDo',
            'elTrigUp', 'elTrigDo',
            'muTrigUp', 'muTrigDo',
            'bTagUp', 'bTagDo',
            'xsecUp', 'xsecDo'
            ]
def getAllVariations() :
    return fakeSystVariations() + mcObjectVariations() + mcWeightVariations()

def fetchVariationHistos(input_fake_file=None, nominal_histo=None, variations=fakeSystVariations()) :
    nom_hname = nominal_histo.GetName()
    return dict([(v, input_fake_file.Get(nom_hname.replace('_NONE','_'+v))) for v in variations])
def computeFakeSysErr2(nominal_histo=None, vars_histos={}) :
    "Compute the bin-by-bin sum2 err including up&down fake systematic variations + stat. unc."
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
def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}
def fetchFakeSysHistosAndComputeSysErr2(input_fake_file=None, nominal_histo=None) :
    vars_histos = fetchVariationHistos(input_fake_file, nominal_histo)
    return computeFakeSysErr2(nominal_histo, vars_histos)
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

