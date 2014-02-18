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
from rootUtils import (buildRatioHistogram,
                       getBinIndices,
                       getNumDenHistos,
                       importRoot)
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptStat(0)
from utils import mkdirIfNeeded

usage="""
Example usage:
%prog \\
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
    plotdir         = opts.plot
    verbose      = opts.verbose
    if verbose : print ('\nUsing the following options:\n'
                        +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))
    fileData = r.TFile.Open(fnameInputDa)
    fileMc   = r.TFile.Open(fnameInputMc)
    if plotdir : mkdirIfNeeded(plotdir)
    assert fileData and fileMc, "Missing input files: data %s, mc %s"%(str(fileData), str(fileMc))
    correctionHistos = {}
    for lep in ['muon', 'elec'] :
        if verbose : print "Lepton: %s"%lep
        hRealDataCr = getNumDenHistos(fileData, lep+'_realCR_all_l_pt')
        hFakeDataLo = getNumDenHistos(fileData, lep+'_fakeHF_all_l_pt')
        hFakeDataHi = getNumDenHistos(fileData, lep+'_fakeHF_high_all_l_pt')
        hFakeMcLo   = getNumDenHistos(fileMc,   lep+'_fakeHF_all_l_pt')
        hFakeMcHi   = getNumDenHistos(fileMc,   lep+'_fakeHF_high_all_l_pt')
        if plotdir :
            hNumDen = [hFakeDataLo, hFakeDataHi, hFakeMcLo, hFakeMcHi]
            for nd in ['num','den'] : plotHistos([h[nd] for h in hNumDen], 'c_'+lep+'_'+nd, plotdir)
            plotHistosRatio(hNumDen, 'c_'+lep+'_ratio', plotdir)
        h2dRealDataCr = getNumDenHistos(fileData, lep+'_realCR_all_l_pt_eta')
        h2dFakeDataLo = getNumDenHistos(fileData, lep+'_fakeHF_all_l_pt_eta')
        h2dFakeDataHi = getNumDenHistos(fileData, lep+'_fakeHF_high_all_l_pt_eta')
        h2dFakeMcLo   = getNumDenHistos(fileMc,   lep+'_fakeHF_all_l_pt_eta')
        h2dFakeMcHi   = getNumDenHistos(fileMc,   lep+'_fakeHF_high_all_l_pt_eta')


        def missingInputHisto(ndHistos) : return any(not h for h in ndHistos.values())
        histoCollToBeChecked = ['hRealDataCr','hFakeDataLo','hFakeDataHi','hFakeMcLo','hFakeMcHi']
        missingHistos = dict([(nhc,hp) for nhc,hp in [(hc, eval(hc)) for hc in histoCollToBeChecked]
                              if missingInputHisto(hp)])
        for v in histoCollToBeChecked :
            print "entries 1d %s : num %d den %d (%s)"%(v, eval(v)['num'].GetEntries(), eval(v)['den'].GetEntries(), str(eval(v)['den']))

        histoCollToBeChecked = ['h2dRealDataCr','h2dFakeDataLo','h2dFakeDataHi','h2dFakeMcLo','h2dFakeMcHi']
        missingHistos = dict([(nhc,hp) for nhc,hp in [(hc, eval(hc)) for hc in histoCollToBeChecked]
                              if missingInputHisto(hp)])

        for v in histoCollToBeChecked :
            print "entries 2d %s : num %d den %d (%s)"%(v, eval(v)['num'].GetEntries(), eval(v)['den'].GetEntries(), str( eval(v)['den']))
        print histoCollToBeChecked
        print missingHistos

        if len(missingHistos) :
            print (lep+' : missing histograms: \n'
                   +'\n'.join(["%s: num %s den %s"%(k, v['num'], v['den'])
                               for k,v in missingHistos.iteritems()]))
            continue


        correctionHistos[lep] = buildCorrectionHisto(hRealDataCr,
                                                     hFakeDataLo, hFakeDataHi,
                                                     hFakeMcLo, hFakeMcHi,
                                                     nIter=nIter, verbose=verbose,
                                                     histoname=lep+'_corHFRate')

        # here do the 2d ones
        print 10*"--"," now doing the 2d ones ",10*"--"
        dummy = h2dRealDataCr['num']
        xAx, yAx = dummy.GetXaxis(), dummy.GetYaxis()
        print dummy.GetName(),": bins (%d, %d)"%(dummy.GetNbinsX(), dummy.GetNbinsY())
        nEtaBins = yAx.GetNbins()
        print 'nEtaBins: ',nEtaBins
        xMin, xMax = xAx.GetXmin(), xAx.GetXmax()
        etaBins = range(1, 1+nEtaBins)
        for eb in etaBins :
            def etaSlice(h, b, p) : return h.ProjectionX(p+h.GetName()+"_eta%d"%b, b, b) # prefix needed to avoid overwriting
            hRealDataCr = dict((k, etaSlice(h, eb, 'rdc')) for k,h in h2dRealDataCr.iteritems())
            hFakeDataLo = dict((k, etaSlice(h, eb, 'fdl')) for k,h in h2dFakeDataLo.iteritems())
            hFakeDataHi = dict((k, etaSlice(h, eb, 'fdh')) for k,h in h2dFakeDataHi.iteritems())
            hFakeMcLo   = dict((k, etaSlice(h, eb, 'fml')) for k,h in h2dFakeMcLo.iteritems())
            hFakeMcHi   = dict((k, etaSlice(h, eb, 'fmh')) for k,h in h2dFakeMcHi.iteritems())
            print "eta bin ",eb
            for k,h in hFakeDataLo.iteritems() :
                print "fakeDataLo %s : %s"%(k, lf2s(binContents(h)))

            correctionHistos[lep+"_eta%d"%eb] = buildCorrectionHisto(hRealDataCr,
                                                                     hFakeDataLo, hFakeDataHi,
                                                                     hFakeMcLo, hFakeMcHi,
                                                                     nIter=nIter, verbose=verbose,
                                                                     histoname=lep+'_corHFRate'+"_eta_bin%d"%eb)
        correctionHistos[lep+"_eta"] = combineEtaSlices(template2d=h2dRealDataCr['num'],
                                                        etaSlicedRates=dict((k,h) for k,h in correctionHistos.iteritems()
                                                                            if (lep+'_eta') in k),
                                                        histoname=lep+'_corHFRate_eta')

        print 10*"--","    done               ",10*"--"
    if verbose : print "saving output to ",fnameOutput
    fileOut = r.TFile.Open(fnameOutput, 'recreate')
    fileOut.cd()
    print 'keys ',correctionHistos.keys()
    for l,h in correctionHistos.iteritems() :
        if verbose : print "%s : writing %s\n%s"%(l, h.GetName(),histo1dToTxt(h))
        h.Write()
    fileOut.Close()



# here be dragons
def buildCorrectionHisto(hndRealDataCr, hndFakeDataLo, hndFakeDataHi, hndFakeMcLo, hndFakeMcHi,
                         histoname='lep_corHFRate', nIter=1, verbose=False) :
    hRealEff = buildRatioHistogram(hndRealDataCr['num'], hndRealDataCr['den'], 'real_eff')
    corrected = dict([(nd, hndFakeDataLo[nd].Clone('corrected_'+nd)) for nd in ['num', 'den']])
    print "buildCorrectionHisto with nIter ",nIter
    can = r.TCanvas('c_'+histoname, '')
    can.Draw()
    can.cd()
    hRealEff.GetYaxis().SetRangeUser(0.0, 1.0)
    hRealEff.Draw('axis')
    tex = r.TLatex()
    tex.SetNDC(True)
    tex.DrawLatex(0.15, 0.915, histoname)
    can._tex = tex
    can._histos = [hRealEff]
    can._leg = r.TLegend(0.925, 0.25, 1.0, 0.9, "iter")
    can._leg.SetBorderSize(0)
    can._leg.SetFillColor(0)
    for iteration in range(nIter) :
        print 'iter ',iteration
        rate = buildRatioHistogram(corrected['num'], corrected['den']) # temporary rate (?)
        if verbose :
            print "Iteration %d, corrected values:"%iteration
            print "  num   %s"%lf2s(binContents(corrected['num']))
            print "  den   %s"%lf2s(binContents(corrected['den']))
            print "  ratio %s"%lf2s(binContents(rate))
            dataNum, dataDen = hndFakeDataHi['num'], hndFakeDataHi['den']
        for nd,tl in [('num','tight'), ('den','loose')] :
            corr, dataLow = corrected[nd], hndFakeDataLo[nd]
            mcLow, mcHi = hndFakeMcLo[nd], hndFakeMcHi[nd]
            corrFact = getCorrFactors(hRealEff, rate, dataNum, dataDen, mcHi, tl)
            corr = correctRate(corr, dataLow, mcLow, corrFact)
        rate.SetLineColor(20+iteration)
        rate.SetMarkerColor(20+iteration)
        rate.SetLineWidth(4)
        rate.SetMarkerStyle(r.kFullCross)
        can._histos.append(rate.DrawClone('el same'))
        can._leg.AddEntry(can._histos[-1], "%d"%iteration, 'lp')
    ratio = buildRatioHistogram(corrected['num'], corrected['den'], histoname)
    can._leg.Draw()
    can.Update()
    can.SaveAs(histoname+'_iterations.png')
    return ratio

def assertSameNbins(histos=[]) :
    nbins = dict([(h.GetName(), h.GetNbinsX()) for h in histos])
    assert len(set(nbins.values()))==1, "different nbin: \n%s"%'\n'.join(["%s : %d"%(kv) for kv in nbin.iteritems()])
def lf2s(l) : return ', '.join(["%.3f"%e for e in l])
def bc(h, b) : return h.GetBinContent(b)
def binContents(h) : return [h.GetBinContent(b) for b in getBinIndices(h)]
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
    bins = getBinIndices(h)
    hisName = h.GetName()
    binEdge = [h.GetBinLowEdge(b) for b in bins]
    binEdge.append(binEdge[-1] + h.GetBinWidth(bins[-1]))
    binCont = binContents(h)
    binErr  = [be(h, b) for b in bins]
    def lf2s(l) : return ', '.join(["%.3f"%e for e in l])

    return '\n'.join(["%s : %s"%(n,v)
                      for n,v in [('hisName',hisName)] + [(l, lf2s(eval(l)))
                                                          for l in ['binEdge', 'binCont', 'binErr']]])

def plotHistos(histos=[], canvasName='c1', outdir='./') :
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
    print "saving ",outdir+'/'+canvasName
    c.SaveAs(outdir+'/'+canvasName+'.png')

def plotHistosRatio(histosPairs=[], canvasName='', outdir='./') :
    histosRatio = [buildRatioHistogram(hh['num'], hh['den']) for hh in histosPairs]
    print histosRatio
    for h in histosRatio : h.SetTitle('ratio '+h.GetTitle())
    plotHistos(histosRatio, canvasName, outdir)
def combineEtaSlices(template2d, etaSlicedRates={}, histoname=''):
    res = template2d.Clone(histoname)
    xAx, yAx = res.GetXaxis(), res.GetYaxis()
    nPtBins, nEtaBins = xAx.GetNbins(), yAx.GetNbins()
    assert nEtaBins==len(etaSlicedRates),"%d eta bins, %d slice histos"%(nEtaBins, len(etaSlicedRates))
    etaBins = range(1, 1+nEtaBins)
    etaHists = [h for k,h in sorted(etaSlicedRates.items())]
    for ie, he in zip(etaBins, etaHists) :
        for b in range(1, 1+nPtBins) :
            res.SetBinContent(b, ie, he.GetBinContent(b))
            res.SetBinError(b, ie, he.GetBinError(b))
    return res


if __name__=='__main__' :
    main()
