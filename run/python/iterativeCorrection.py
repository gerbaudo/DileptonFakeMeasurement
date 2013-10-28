#!/bin/env python

# Determine the contamination from prompt leptons to the HF estimate
# in the high-pt tail.

# The procedure is based on the idea described in section 7.1.3 of
# ATL-COM-PHYS-2011-768.
# This python implementation is based on the c++ implementation by
# Matt (mrelich6@gmail.com), originally in IterativeFakeCorrection.cxx
#
# Here is a brief description of the idea:
# - For the matrix method we need to compute p(T|R) and p(T|F); we
#   compute them as N_T/N_L in each one of two cases (R or F)
#   as a function of p_T
# - N_T and N_L are taken from data. For the HF fakes, we do a
#   tag-and-probe in a ccbar/bbar region (see
#   MeasureFakeRate2::passHFCR)
# - However, N_T and N_L are contaminated by real, prompt leptons from
#   ewk production. So what we do is to subtract the real
#   contribution. This is done by identifying two regions in m_T:
#   CR_HF with m_T<40 (mostly fake), and CR_HF_high with m_T<100
#   (where there is more contamination from real leptons).
# - 0^th iteration: p(T|F) = (  N_T(CR_HF)
#                             / N_L(CR_HF))
# - 1^st iteration: p(T|F) = (  (N_T(CR_HF) - c^1_T * N_T(CR_HF, MC))
#                             / (N_L(CR_HF) - c^1_L * N_L(CR_HF, MC)))
#   where c^1_T = (  (data(CR_HF_high) - fake(pred, CR_HF_high))
#                  / mc(pred. all-mc-except-hf, CR_HF_high))
# - at the following iterations, you update the fake estimate in the
#   numerator of the correction factor c with the fake prediction that
#   you can now compute with the one-lepton matrix method.

# davide.gerbaudo@gmail.com
# October 2013

from math import sqrt
import optparse
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options

usage="""
Example usage:
IterFake \\
 --input_data out/fakerate/merged/data_${TAG}.root \\
 --input_mc   out/fakerate/merged/allBkgButHf_${TAG}.root \\
 --output     out/fakerate/merged/iterative_out_${TAG}.root \\
 >& log/fakerate/IterFake_${TAG}.log
"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-n', '--n_iter', type='int', default=8)
    parser.add_option('-m', '--input_mc')
    parser.add_option('-d', '--input_data')
    parser.add_option('-o', '--output')
    parser.add_option('-p', '--plot', help='plot inputs') # todo: implement sanity plot vs. n_iter
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['n_iter', 'input_mc', 'input_data', 'output']
    otherOptions = ['plot', 'verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    nIter        = opts.n_iter
    fnameInputMc = opts.input_mc
    fnameInputDa = opts.input_data
    fnameOutput  = opts.output
    plot         = opts.plot
    verbose      = opts.verbose
    if verbose : print ('\nUsing the following options:\n'
                        +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))
    fileData = r.TFile.Open(fnameInputDa)
    fileMc   = r.TFile.Open(fnameInputMc)
    assert fileData and fileMc, "Missing input files: data %s, mc %s"%(str(fileData), str(fileMc))
    def getNumDenHistos(f, baseHname='base_histo_name', suffNum='_num', suffDen='_den') :
        "baseHname is something like 'lep_controlreg_chan_var', see MeasureFakeRate2::initHistos()"
        num = f.Get(baseHname+suffNum)
        den = f.Get(baseHname+suffDen)
        return {'num':num, 'den':den}
    correctionHistos = {}
    for lep in ['muon', 'elec'] :
        if verbose : print "Lepton: %s"%lep
        hRealDataCr = getNumDenHistos(fileData, lep+'_realCR_all_l_pt')
        hFakeDataLo = getNumDenHistos(fileData, lep+'_fakeHF_all_l_pt')
        hFakeDataHi = getNumDenHistos(fileData, lep+'_fakeHF_high_all_l_pt')
        hFakeMcLo   = getNumDenHistos(fileMc,   lep+'_fakeHF_all_l_pt')
        hFakeMcHi   = getNumDenHistos(fileMc,   lep+'_fakeHF_high_all_l_pt')
        if plot :
            hNumDen = [hFakeDataLo, hFakeDataHi, hFakeMcLo, hFakeMcHi]
            for nd in ['num','den'] : plotHistos([h[nd] for h in hNumDen], 'c_'+lep+'_'+nd)
            plotHistosRatio(hNumDen, 'c_'+lep+'_ratio')
        def missingInputHisto(ndHistos) : return any(not h for h in ndHistos.values())
        histoCollToBeChecked = ['hRealDataCr','hFakeDataLo','hFakeDataHi','hFakeMcLo','hFakeMcHi']
        missingHistos = dict([(nhc,hp) for nhc,hp in [(hc, eval(hc)) for hc in histoCollToBeChecked]
                              if missingInputHisto(hp)])
        if len(missingHistos) :
            print (lep+' : missing histograms: \n'
                   +'\n'.join(["%s: num %s den %s"%(k, v['num'], v['den'])
                               for k,v in missingHistos.iteritems()]))
            continue
        hRealEff = ratioHistogram(hRealDataCr['num'], hRealDataCr['den'], 'real_eff')
        corrected = dict([(nd, hFakeDataLo[nd].Clone('corrected_'+nd)) for nd in ['num', 'den']])
        for iteration in range(nIter) :
            rate = ratioHistogram(corrected['num'], corrected['den']) # temporary rate (?)
            if verbose :
                def lf2s(l) : return ', '.join(["%.3f"%e for e in l])
                print "Iteration %d, corrected values:"%iteration
                print "  num   %s"%lf2s(binContents(corrected['num']))
                print "  den   %s"%lf2s(binContents(corrected['den']))
                print "  ratio %s"%lf2s(binContents(rate))
            dataNum, dataDen = hFakeDataHi['num'], hFakeDataHi['den']
            for nd,tl in [('num','tight'), ('den','loose')] :
                corr, dataLow = corrected[nd], hFakeDataLo[nd]
                mcLow, mcHi = hFakeMcLo[nd], hFakeMcHi[nd]
                corrFact = getCorrFactors(hRealEff, rate, dataNum, dataDen, mcHi, tl)
                corr = correctRate(corr, dataLow, mcLow, corrFact)
        ratio = ratioHistogram(corrected['num'], corrected['den'], lep+'_corHFRate')
        correctionHistos[lep] = ratio
    if verbose : print "saving output to ",fnameOutput
    fileOut = r.TFile.Open(fnameOutput, 'recreate')
    fileOut.cd()
    for l,h in correctionHistos.iteritems() :
        if verbose : print "%s : writing %s\n%s"%(l, h.GetName(),histo1dToTxt(h))
        h.Write()
    fileOut.Close()

def assertSameNbins(histos=[]) :
    nbins = dict([(h.GetName(), h.GetNbinsX()) for h in histos])
    assert len(set(nbins.values()))==1, "different nbin: \n%s"%'\n'.join(["%s : %d"%(kv) for kv in nbin.iteritems()])
def bc(h, b) : return h.GetBinContent(b)
def binContents(h) : return [h.GetBinContent(b) for b in range(1, 1+h.GetNbinsX())]
def be(h, b) : return h.GetBinError(b)
def getCorrFactors(hRealEff, hRate, hDataNum, hDataDen, hMc, tightOrLoose) :
    """
    Bin-by-bin correction factor:
    C = (N^(data,high) - N^(fake pred, high))/N^(MC,high)
    """
    assertSameNbins([hRealEff, hRate, hDataNum, hDataDen, hMc])
    isTight = tightOrLoose=='tight'
    def getFakePred(real, fake, loose, tight) : return ( (fake if isTight else 1.0)
                                                         / (real-fake)
                                                         * (loose*real - tight) )
    bins = range(1, 1+hRate.GetNbinsX())
    loose = [bc(hDataDen, b) for b in bins]
    tight = [bc(hDataNum, b) for b in bins]
    data  = tight if isTight else loose
    mc    = [bc(hMc,      b) for b in bins]
    real  = [bc(hRealEff, b) for b in bins]
    fake  = [bc(hRate,    b) for b in bins]
    fakeP = [getFakePred(r, f, l, t) for r, f, l, t in zip(real, fake, loose, tight)]
    return [(d - f)/m for d,f,m in zip(data, fakeP, mc)]
def correctRate(hRate, hData, hMc, corrections) :
    assertSameNbins([hRate, hData, hMc])
    bins = range(1, 1+hRate.GetNbinsX())
    assert len(bins)==len(corrections)," bins[%d], corr[%d]"%(len(bins), len(corrections))
    daCnt, daErr = [bc(hData, b) for b in bins], [be(hData, b) for b in bins]
    mcCnt, mcErr = [bc(hMc, b)   for b in bins], [be(hMc, b)   for b in bins]
    def sumquad(a, b) : return sqrt(a*a + b*b)
    for b, c, dc, de, mc, me in zip(bins, corrections, daCnt, daErr, mcCnt, mcErr) :
        hRate.SetBinContent(b, dc - c * mc)
        hRate.SetBinError  (b, sumquad(de, c*me))
    return hRate
def histo1dToTxt(h) :
    "represent a TH1 as a string with name, bin edges, contents, and errors"
    bins = range(1, 1+h.GetNbinsX())
    hisName = h.GetName()
    binEdge = [h.GetBinLowEdge(b) for b in bins]
    binEdge.append(binEdge[-1] + h.GetBinWidth(bins[-1]))
    binCont = binContents(h)
    binErr  = [be(h, b) for b in bins]
    def lf2s(l) : return ', '.join(["%.3f"%e for e in l])
    
    return '\n'.join(["%s : %s"%(n,v)
                      for n,v in [('hisName',hisName)] + [(l, lf2s(eval(l)))
                                                          for l in ['binEdge', 'binCont', 'binErr']]])
def ratioHistogram(num, den, name='ratio') :
    r = num.Clone(name)
    r.SetDirectory(0) # we usually don't care about the ownership of these temporary objects
    r.Reset()
    r.Divide(num, den, 1, 1)
    return r

def plotHistos(histos=[], canvasName='c1') :
    c = r.TCanvas(canvasName)
    c.cd()
    maximum = max([h.GetMaximum() for h in histos])
    padMaster = None
    colors = [r.kBlack, r.kRed, r.kBlue, r.kGreen, r.kViolet]
    markers = [r.kOpenCircle, r.kOpenSquare, r.kOpenTriangleUp, r.kOpenDiamond, r.kOpenCross]
    for h,col,m in zip(histos, colors, markers) :
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerStyle(m)
        h.Draw('same' if padMaster else '')
        padMaster = padMaster if padMaster else h
    padMaster.SetMaximum(1.1*maximum)
    leg = c.BuildLegend()
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    c.Update()
    c.SaveAs(canvasName+'.png')

def plotHistosRatio(histosPairs=[], canvasName='') :
    histosRatio = [ratioHistogram(hh['num'], hh['den']) for hh in histosPairs]
    plotHistos(histosRatio, canvasName)

if __name__=='__main__' :
    main()
