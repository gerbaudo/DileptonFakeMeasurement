#!/bin/env python

# Determine the contamination from prompt leptons to the HF estimate
# in the high-pt tail.

# The procedure is based on the idea described in section 7.1.3 of
# ATL-COM-PHYS-2011-768.
# This python implementation is based on the c++ implementation by
# Matt (mrelich6@gmail.com), originally in IterativeFakeCorrection.cxx

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
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['n_iter', 'input_mc', 'input_data', 'output']
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    nIter        = opts.n_iter
    fnameInputMc = opts.input_mc
    fnameInputDa = opts.input_data
    fnameOutput  = opts.output
    verbose      = opts.verbose
    if verbose : print ('\nUsing the following options:\n'
                        +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in requiredOptions))
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
        hFakeDataLo = getNumDenHistos(fileData, lep+'_fakeHF_all_l_pt') # bug? shouldnt be 'low'?
        hFakeDataHi = getNumDenHistos(fileData, lep+'_fakeHF_high_all_l_pt')
        hFakeMcLo   = getNumDenHistos(fileMc,   lep+'_fakeHF_all_l_pt') # bug? shouldnt be 'low'?
        hFakeMcHi   = getNumDenHistos(fileMc,   lep+'_fakeHF_high_all_l_pt')        
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
    if verbose :
        print "saving output to ",fnameOutput
    fileOut = r.TFile.Open(fnameOutput, 'recreate')
    fileOut.cd()
    for l,h in correctionHistos.iteritems() :
        if verbose :
            print "%s : writing %s"%(l,h.GetName())
            print histo1dToTxt(h)
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

if __name__=='__main__' :
    main()
