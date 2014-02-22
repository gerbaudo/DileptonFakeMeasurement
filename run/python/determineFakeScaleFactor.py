#!/bin/env python

# Determine the scale factors needed for the fake estimate.

# inputs:
# 'all_l_pt' histograms for a few control regions, (fill in details here)
# from data_${TAG}.root, allBkg_${TAG}.root
#
# outputs:
# 5 scale float factors : mu_[real, hf], el[real, hf, conv]

import datetime
from math import sqrt
import optparse
from pdgRounding import pdgRound
from rootUtils import (getNumDenHistos
                       ,buildRatioHistogram
                       ,buildBotTopPads
                       ,drawLegendWithDictKeys
                       ,importRoot
                       )
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)
from utils import (first
                   ,rmIfExists
                   ,mkdirIfNeeded
                   )

usage="""
Example usage:
%prog \\
 --tag ${TAG} \\
 --input_dir  out/fakerate/merged/ \\
 --input_iter out/fakerate/merged/iterative_out_${TAG}.root \\
 --output_dir out/fakerate/merged/plot_${TAG} \\
 >& log/fakerate/FakePlot_${TAG}.log
"""
def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-c', '--input_iter')
    parser.add_option('-o', '--output_dir')
    parser.add_option('-O', '--output_file',help='store ratio histograms here')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_dir', 'input_iter', 'output_dir']
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag            = opts.tag.strip('_')
    inputDirname   = opts.input_dir
    fnameInputIter = opts.input_iter
    outputDirname  = opts.output_dir
    outputDirname  = outputDirname+'/' if not outputDirname.endswith('/') else outputDirname
    outputFilename = opts.output_file
    mkdirIfNeeded(outputDirname)
    verbose        = opts.verbose
    if verbose : print ('\nUsing the following options:\n'
                        +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))
    inputDirname = inputDirname+'/' if not inputDirname.endswith('/') else inputDirname
    fileData = r.TFile.Open(inputDirname+'data_'       +tag+'.root')
    fileMc   = r.TFile.Open(inputDirname+'allBkg_'     +tag+'.root')
    fileHf   = r.TFile.Open(inputDirname+'heavyflavor_'+tag+'.root')
    fileIter = r.TFile.Open(fnameInputIter)
    assert fileData, "Missing input file data %s"%str(fileData)
    assert fileMc,   "Missing input file mc   %s"%str(fileMc)
    assert fileHf,   "Missing input file hf   %s"%str(fileHf)
    assert fileIter, "Missing input file iter %s"%str(fileIter)
    outputFile = r.TFile.Open(outputFilename, 'recreate') if outputFilename else None
    el_conv_sf = computeAndPlotConvSf(fileData, fileMc, 'elec', 'all_l_pt', outputDirname, outputFile)
    el_qcd_sf  = computeAndPlotHfSf  (fileIter, fileHf, 'elec', 'all_l_pt', outputDirname, outputFile)
    mu_qcd_sf  = computeAndPlotHfSf  (fileIter, fileHf, 'muon', 'all_l_pt', outputDirname, outputFile)
    el_real_sf = computeAndPlotRealSf(fileData, fileMc, 'elec', 'all_l_pt', outputDirname)
    mu_real_sf = computeAndPlotRealSf(fileData, fileMc, 'muon', 'all_l_pt', outputDirname)

    el_conv_sf2d = computeAndPlotConvSf2d(fileData, fileMc, 'elec', 'all_l_pt', outputDirname)
    el_qcd_sf2d  = computeAndPlotHfSf2d  (fileIter, fileHf, 'elec', 'all_l_pt', outputDirname)
    mu_qcd_sf2d  = computeAndPlotHfSf2d  (fileIter, fileHf, 'muon', 'all_l_pt', outputDirname)

    print "# --- paste the lines below in buildWeightedMatrix.py ---"
    print "# %s, %s"%(tag, datetime.datetime.now())
    print "mu_qcdSF, mu_realSF = %s, %s"%(mu_qcd_sf, mu_real_sf)
    print "el_convSF, el_qcdSF, el_realSF = %s, %s, %s"%(el_conv_sf, el_qcd_sf, el_real_sf)
    if outputFile : outputFile.Close()

def computeAndPlotConvSf(fileData, fileMc, lepton, variable_name, outdir, outfile=None) :
    "Electron conversion: simplest case, just data/mc"
    eff_da = buildRate(fileData, lepton+'_fakeConv_'+variable_name)
    eff_mc = buildRate(fileMc,   lepton+'_fakeConv_'+variable_name)
    ratio  = buildRatioHistogram(eff_da, eff_mc)
    fitFunc = fitWithConst(ratio)
    p0, p0Err, chi2, ndf = fitResults(fitFunc)
    p0, p0Err = pdgRound(p0, p0Err)
    print "SF for %s fake conv : %s +/- %s"%(lepton, p0, p0Err)
    graphics = {'xtitle' : xTitle(lepton, variable_name),
                'ytitle' : lepton+' p(tight | fake conv)',
                'colors' : {'data' : r.kBlack, 'mc' : mcColor(lepton)},
                'markers': {'data' : r.kFullCircle, 'mc' : mcMarker(lepton)},
                'labels' : {'data' : 'Data: Conversion CR',
                            'mc'   : 'MC Comb: Conv CR'}}
    plotHistRatioAndFit({'data':eff_da, 'mc':eff_mc}, ratio, fitFunc, outdir+'/fit_'+lepton+'_conv', graphics)
    if outfile : saveObject(outfile, ratio, 'elec_convSF_pt')
    return p0

def computeAndPlotConvSf2d(fileData, fileMc, lepton, variable_name, outdir) :
    "Electron conversion: simplest case, just data/mc"
    histoname = 'elec_fakeConv_all_l_pt_eta'
    eff_da = buildRate(fileData, histoname)
    eff_mc = buildRate(fileMc,   histoname)
    print 'SF conversion: using histo ',histoname
    ratio  = buildRatioHistogram(eff_da, eff_mc)
    ratio.Print()
    xAx, yAx = ratio.GetXaxis(), ratio.GetYaxis()
    #pt_eta : check that x is pt, y is eta
    nEtaBins = yAx.GetNbins()
    xMin, xMax = xAx.GetXmin(), xAx.GetXmax()
    fitFunc = r.TF1('fit_func_const_'+ratio.GetName(), '[0]', xMin, xMax)
    etaBins = range(1, 1+nEtaBins)
    slices = [ratio.ProjectionX("%s_bin%d"%(ratio.GetName(), b), b, b, 'e') for b in etaBins]
    for b, s in zip(etaBins, slices) :
        s.SetTitle("data/mc conversion : eta bin %d"%b)
        s.Fit(fitFunc.GetName(), '0RQ') # do not draw, range, quiet
        p0, p0Err, chi2, ndf = fitResults(fitFunc)
        p0, p0Err = pdgRound(p0, p0Err)
        print "bin %d :  %s +/- %s"%(b, p0, p0Err)
        can = r.TCanvas('')
        s.Draw('ep')
        fitFunc.Draw('same')
        tex = r.TLatex()
        tex.SetNDC(True)
        fitParLabel = "Const. fit : %s #pm %s"%(p0, p0Err)
        fitGoodLabel = "#chi^{2}/DOF : %.2f / %d"%(chi2, ndf)
        tex.SetTextSize(yAx.GetTitleSize())
        tex.SetTextFont(yAx.GetTitleFont())
        tex.DrawLatex(0.15, 0.45, s.GetTitle())
        tex.DrawLatex(0.15, 0.40, "#splitline{%s}{%s}"%(fitParLabel, fitGoodLabel))
        can.Update()
        for ext in ['eps','png'] : can.SaveAs(outdir+"/fit_el_conv_etabin%d.%s"%(b, ext))
#     graphics = {'xtitle' : xTitle(lepton, variable_name),
#                 'ytitle' : lepton+' p(tight | fake conv)',
#                 'colors' : {'data' : r.kBlack, 'mc' : mcColor(lepton)},
#                 'markers': {'data' : r.kFullCircle, 'mc' : mcMarker(lepton)},
#                 'labels' : {'data' : 'Data: Conversion CR',
#                             'mc'   : 'MC Comb: Conv CR'}}
#     plotHistRatioAndFit({'data':eff_da, 'mc':eff_mc}, ratio, fitFunc, outdir+lepton+'_fakeconv',
#                         graphics)
    return p0

def computeAndPlotHfSf(fileIter, fileHf, lepton, variable_name, outdir, outfile=None) :
    "HF tag and probe; in this case we need to subract out the contamination"
    eff_da = fileIter.Get(lepton+'_corHFRate')
    eff_mc = buildRate(fileHf, lepton+'_fakeHF_'+variable_name)
    ratio = buildRatioHistogram(eff_da, eff_mc)
    fitFunc = fitWithConst(ratio)
    p0, p0Err, chi2, ndf = fitResults(fitFunc)
    p0, p0Err = pdgRound(p0, p0Err)
    print "SF for %s fake HF : %s +/- %s"%(lepton, p0, p0Err)
    graphics = {'xtitle' : xTitle(lepton, variable_name),
                'ytitle' : lepton+' p(tight | fake hf)',
                'colors' : {'data' : r.kBlack, 'mc' : mcColor(lepton)},
                'markers': {'data' : r.kFullCircle, 'mc' : mcMarker(lepton)},
                'labels' : {'data' : 'Data HF Tag and Probe (Iterative Subtraction)',
                            'mc'   : 'b#bar{b}/c#bar{c} MC: HF Tag and Probe'}}
    plotHistRatioAndFit({'data':eff_da, 'mc':eff_mc}, ratio, fitFunc, outdir+'/fit_'+lepton+'_fakehf',
                        graphics)
    if outfile : saveObject(outfile, ratio, lepton+'_qcdSF_pt')
    return p0
def computeAndPlotHfSf2d(fileIter, fileHf, lepton, variable_name, outdir) :
    eff_da = fileIter.Get(lepton+'_corHFRate_eta')
    eff_mc = buildRate(fileHf, lepton+'_fakeHF_all_l_pt_eta')
    ratio = buildRatioHistogram(eff_da, eff_mc)
    ratio.Print()
    xAx, yAx = ratio.GetXaxis(), ratio.GetYaxis()
    print ratio.GetName(),": bins (%d, %d)"%(ratio.GetNbinsX(), ratio.GetNbinsY())
    nEtaBins = yAx.GetNbins()
    print 'nEtaBins: ',nEtaBins
    xMin, xMax = xAx.GetXmin(), xAx.GetXmax()
    fitFunc = r.TF1('fit_func_const_'+ratio.GetName(), '[0]', xMin, xMax)
    etaBins = range(1, 1+nEtaBins)
    slices = [ratio.ProjectionX("%s_bin%d"%(ratio.GetName(), b), b, b, 'e') for b in etaBins]
    for b, s in zip(etaBins, slices) :
        s.SetTitle(lepton+" data/mc heavyflavor : eta bin %d"%b)
        s.Fit(fitFunc.GetName(), '0RQ') # do not draw, range, quiet
        p0, p0Err, chi2, ndf = fitResults(fitFunc)
        p0, p0Err = pdgRound(p0, p0Err)
        print "bin %d :  %s +/- %s"%(b, p0, p0Err)
        can = r.TCanvas('')
        s.Draw('ep')
        fitFunc.Draw('same')
        tex = r.TLatex()
        tex.SetNDC(True)
        fitParLabel = "Const. fit : %s #pm %s"%(p0, p0Err)
        fitGoodLabel = "#chi^{2}/DOF : %.2f / %d"%(chi2, ndf)
        tex.SetTextSize(yAx.GetTitleSize())
        tex.SetTextFont(yAx.GetTitleFont())
        tex.DrawLatex(0.15, 0.45, s.GetTitle())
        tex.DrawLatex(0.15, 0.40, "#splitline{%s}{%s}"%(fitParLabel, fitGoodLabel))
        can.Update()
        for ext in ['eps','png'] : can.SaveAs(outdir+'/fit_'+lepton+"_heavyflavor_etabin%d.%s"%(b, ext))
    return p0

def computeAndPlotRealSf(file_data, file_mc, lepton, variable_name, outdir) :
    "Scale factor from the real control region, Z tag and probe"
    eff_da = buildSideBandSubRate(file_data, lepton, variable_name)
    eff_mc = buildSideBandSubRate(file_mc,   lepton, variable_name)
    ratio = buildRatioHistogram(eff_da, eff_mc)
    fitFunc = fitWithConst(ratio)
    p0, p0Err, chi2, ndf = fitResults(fitFunc)
    p0, p0Err = pdgRound(p0, p0Err)
    print "SF for %s real : %s +/- %s"%(lepton, p0, p0Err)
    graphics = {'xtitle' : xTitle(lepton, variable_name),
                'ytitle' : lepton+' p(tight | real)',
                'colors' : {'data' : r.kBlack, 'mc' : mcColor(lepton)},
                'markers': {'data' : r.kFullCircle, 'mc' : mcMarker(lepton)},
                'labels' : {'data' : 'Data: Z Tag and Probe',
                            'mc'   : 'MC Comb: Z Tag and Probe'}}
    plotHistRatioAndFit({'data':eff_da, 'mc':eff_mc}, ratio, fitFunc, outdir+'/fit_'+lepton+'_real',
                        graphics)
    return p0

def buildRate(file, histo_basename) :
    hs = getNumDenHistos(file, histo_basename)
    return buildRatioHistogram(hs['num'], hs['den'])
def fitWithConst(histo) :
    xMin, xMax = histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax()
    fitFunc = r.TF1('fit_func_const_'+histo.GetName(), '[0]', xMin, xMax)
    histo.Fit(fitFunc.GetName(), '0RQ') # do not draw, range, quiet
    return fitFunc
def fitResults(fitFunc) :
    "read the fit parameters from the function"
    p0, p0Err = fitFunc.GetParameter(0), fitFunc.GetParError(0)
    chi2, ndf = fitFunc.GetChisquare(), fitFunc.GetNDF()
    return p0, p0Err, chi2, ndf
def buildSideBandSubRate(file, lepton, variable_name) :
    hZwindow = getNumDenHistos(file, lepton+'_realCR_'      +variable_name)
    hSideLo  = getNumDenHistos(file, lepton+'_realSideLow_' +variable_name)
    hSideHi  = getNumDenHistos(file, lepton+'_realSideHigh_'+variable_name)
    def sbErr(el, ew, eh) :
        "SidebandError; DG this doesn't make any sense to me; ask Matt"
        return sqrt(ew*ew + (el*el + eh*eh))
    def be(h) : return [h.GetBinError(b) for b in range(1, 1+h.GetNbinsX())]
    errs = dict([(k, [sbErr(l,w,h)
                      for l, w, h in zip(be(hSideLo[k]), be(hZwindow[k]), be(hSideHi[k]))])
                 for k in ['num','den']])
    hSideLo ['num'].Add(hSideHi['num'])
    hSideLo ['den'].Add(hSideHi['den'])
    hZwindow['num'].Add(hSideLo['num'], -1)
    hZwindow['den'].Add(hSideLo['den'], -1)
    for nd in ['num', 'den'] :
        h, err = hZwindow[nd], errs[nd]
        nbins = h.GetNbinsX()
        for i, b in  zip(range(nbins), range(1, 1+nbins)) :
            h.SetBinError(b, err[i])
    return buildRatioHistogram(hZwindow['num'], hZwindow['den'])
def plotHistRatioAndFit(histos, ratio, fitfunc, outfname, graphics_attributes={}) :
    can = r.TCanvas('can_'+outfname, outfname, 800, 600)
    botPad, topPad = buildBotTopPads(can, splitFraction=0.35)
    can.cd()
    topPad.Draw()
    topPad.cd()
    pm = first(histos)
    pm.Draw('axis')
    xAx, yAx = pm.GetXaxis(), pm.GetYaxis()
    xAx.SetTitle('')
    xAx.SetLabelSize(0)
    yAx.SetRangeUser(0.0, 1.2)
    yAx.SetNdivisions(505) # for some reason some of the inputs are already ratios with -201
    textScaleUp = 1.0/topPad.GetHNDC()
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetTitle(graphics_attributes['ytitle'])
    yAx.SetTitleOffset(yAx.GetTitleOffset()/textScaleUp)
    for k, h in histos.iteritems() :
        h.SetMarkerStyle(graphics_attributes['markers'][k])
        h.SetMarkerColor(graphics_attributes['colors'][k]) 
        h.SetLineColor(graphics_attributes['colors'][k])
        h.Draw('same')
    labels = graphics_attributes['labels']
    drawLegendWithDictKeys(topPad, dict([(labels[s], histos[s]) for s in histos.keys()]), legWidth=0.4)
    can.cd()
    botPad.Draw()
    botPad.cd()
    ratio.Draw('axis')
    textScaleUp = 1.0/botPad.GetHNDC()
    xAx, yAx = ratio.GetXaxis(), ratio.GetYaxis()
    yAx.SetRangeUser(0.0, 2.0)
    xAx.SetTitle(graphics_attributes['xtitle'])
    yAx.SetNdivisions(-202)
    yAx.SetTitle('Data/Sim')
    yAx.CenterTitle()
    for a in [xAx, yAx] :
        if a.GetLabelSize()>0.04 : continue
        a.SetLabelSize(a.GetLabelSize()*textScaleUp)
        a.SetTitleSize(a.GetTitleSize()*textScaleUp)
        a.SetTitleOffset(a.GetTitleOffset()/textScaleUp)
    ratio.SetMarkerStyle(graphics_attributes['markers']['data'])
    fitfunc.SetLineWidth(2)
    fitfunc.SetLineStyle(2)
    fitfunc.Draw('same')
    ratio.Draw('same')
    tex = r.TLatex()
    tex.SetNDC(True)
    p0, p0Err, chi2, ndf = fitResults(fitfunc)
    fitParLabel = "Const. fit : %s #pm %s"%(pdgRound(p0, p0Err))
    fitGoodLabel = "#chi^{2}/DOF : %.2f / %d"%(chi2, ndf)
    tex.SetTextSize(yAx.GetTitleSize())
    tex.SetTextFont(yAx.GetTitleFont())
    tex.DrawLatex(0.15, 0.40, "#splitline{%s}{%s}"%(fitParLabel, fitGoodLabel))
    can.Update()
    for ext in ['eps','png'] :
        rmIfExists(outfname+'.'+ext) # avoid root warnings
        can.SaveAs(outfname+'.'+ext)
def mcColor(lepton) : return {'elec':r.kRed, 'muon':r.kBlue}[lepton]
def mcMarker(lepton) : return {'elec':r.kOpenSquare, 'muon':r.kOpenTriangleUp}[lepton]
def xTitle(lepton, variable_name) :
    return lepton +' p_{T} [GeV]' if 'l_pt' in variable_name else ''
def saveObject(dest, obj, name) :
    oPwd = r.gDirectory
    dest.cd()
    obj.Write(name)
    oPwd.cd()

if __name__=='__main__' :
    main()
