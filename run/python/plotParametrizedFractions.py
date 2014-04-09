#!/bin/env python

# Plot the fake lepton composition as a function of several variables

# davide.gerbaudo@gmail.com
# Feb 2014


import optparse
from rootUtils import (importRoot
                       ,buildRatioHistogram
                       ,drawLegendWithDictKeys
                       ,getMinMax
                       ,getBinContents
                       ,getBinIndices
                       ,topRightLegend
                       )
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
from fakeUtils import (samples
                       ,fakeProcesses
                       ,getInputFiles
                       ,buildRatio
                       )
from utils import (first
                   ,commonPrefix
                   ,commonSuffix
                   ,mkdirIfNeeded
                   ,rmIfExists
                   )
import SampleUtils

usage="""
Example usage:
%prog \\
 --tag ${TAG} \\
 --input_dir out/fakerate/merged/ \\
 --output_file out/fakerate/merged/fraction_${TAG}.root \\
 --output_plot out/fakerate/merged/fraction_plots_${TAG} \\
 >& log/fakerate/fractions.log
"""
def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-o', '--output_file')
    parser.add_option('-p', '--output_plot')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_dir', 'output_file', 'output_plot']
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag = opts.tag
    inputDirname  = opts.input_dir
    outputFname   = opts.output_file
    outputPlotDir = opts.output_plot
    verbose       = opts.verbose
    if verbose : print '\nUsing the following options:\n'+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

    allInputFiles = getInputFiles(inputDirname, tag, verbose) # includes allBkg, which is used only for sys
    assert all(f for f in allInputFiles.values()), ("missing inputs: \n%s"%'\n'.join(["%s : %s"%kv for kv in allInputFiles.iteritems()]))
    outputPlotDir = outputPlotDir+'/' if not outputPlotDir.endswith('/') else ''
    mkdirIfNeeded(outputPlotDir)
    inputFiles = dict((k, v) for k, v in allInputFiles.iteritems() if k in fakeProcesses())
    if verbose : print 'Using the following input files:\n'+'\n'.join(["%s : %s"%(p, f.GetName()) for p,f in inputFiles.iteritems()])
    outputFile = r.TFile.Open(outputFname, 'recreate')
    buildMuonFractions    (inputFiles, outputPlotDir, outputFile, verbose)
    buildElectronFractions(inputFiles, outputPlotDir, outputFile, verbose)
    outputFile.Close()
    if verbose : print "output plots saved to %s"%outputPlotDir
    if verbose : print "output histos saved to %s"%outputFname

def independentVariables() :
    """These are the variables vs. which we filled the flavor
    histograms; the 2d histo name is something like
        <prefix>_flavor_<var>
    See MeasureFakeRate2::initHistos()
    """
    return ["pt", "pt_etaC", "pt_etaF", "eta", "metrel"]

def regionsToBePlotted() :
    "regions for which we plot the fractions"
    return  ['CR_SRWHnoMlj', 'CR_SRWH1j']
    #return  ['CR_WHZVfake1j', 'CR_WHZVfake2j', 'CR_WHfake1j', 'CR_WHfake2j', 'CR_WHZV1j', 'CR_WHZV2j', 'CR_SRWH1j', 'CR_SRWH2j', 'CR_SRWHnoMlj']


def buildMuonFractions(inputFiles, outplotdir, outfile, verbose=False) :
    ""
    for reg in regionsToBePlotted() :
        for var in independentVariables() :

            binLabel = 'real' # real fractions
            # draw inputs
            frameTitle = 'real muon: '+reg+' loose;'+var
            histoName = 'muon_'+reg+'_all_flavor_'+var+'_den'
            canvasName = 'muon_real_'+reg+'_'+var+'_den'
            histosFlavorSliceReal = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSliceReal, canvasName, outplotdir, frameTitle)

            # normalize and draw fractions (den only)
            normalizeHistos(histosFlavorSliceReal)
            frameTitle = 'real muon: '+reg+';'+var
            canvasBaseName = 'muon_real_'+reg+'_'+var+'_frac'
            plotFractionsHistos (histosFlavorSliceReal, canvasBaseName+'_histo', outplotdir, frameTitle)
            plotFractionsStacked(histosFlavorSliceReal, canvasBaseName+'_stack', outplotdir, frameTitle)
            saveHistos(histosFlavorSliceReal.values(), outfile)

            binLabel = 'qcd'
            # draw inputs, num and den
            frameTitle = 'qcd muon: '+reg+' tight;'+var
            histoName = 'muon_'+reg+'_all_flavor_'+var+'_num'
            canvasName = 'muon_qcd_'+reg+'_'+var+'_num'
            histosFlavorSlice = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSlice, canvasName, outplotdir, frameTitle)

            frameTitle = 'qcd muon: '+reg+' loose;'+var
            histoName = 'muon_'+reg+'_all_flavor_'+var+'_den'
            canvasName = 'muon_qcd_'+reg+'_'+var+'_den'
            histosFlavorSlice = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSlice, canvasName, outplotdir, frameTitle)

            # normalize and draw fractions (den only)
            normalizeHistos(histosFlavorSlice)
            frameTitle = 'qcd muon: '+reg+';'+var
            canvasBaseName = 'muon_qcd_'+reg+'_'+var+'_frac'
            plotFractionsHistos (histosFlavorSlice, canvasBaseName+'_histo', outplotdir, frameTitle)
            plotFractionsStacked(histosFlavorSlice, canvasBaseName+'_stack', outplotdir, frameTitle)
            saveHistos(histosFlavorSlice.values(), outfile)

def buildElectronFractions(inputFiles, outplotdir, outfile, verbose=False) :
    ""
    for reg in regionsToBePlotted() :
        for var in independentVariables() :
            binLabel = 'real' # real fractions
            # draw inputs
            frameTitle = 'real elec: '+reg+' loose;'+var
            histoName = 'elec_'+reg+'_all_flavor_'+var+'_den'
            canvasName = 'elec_real_'+reg+'_'+var+'_den'
            histosFlavorSliceReal = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSliceReal, canvasName, outplotdir, frameTitle)

            # normalize and draw fractions (den only)
            normalizeHistos(histosFlavorSliceReal)
            frameTitle = 'real elec: '+reg+';'+var
            canvasBaseName = 'elec_real_'+reg+'_'+var+'_frac'
            plotFractionsHistos (histosFlavorSliceReal, canvasBaseName+'_histo', outplotdir, frameTitle)
            plotFractionsStacked(histosFlavorSliceReal, canvasBaseName+'_stack', outplotdir, frameTitle)
            saveHistos(histosFlavorSliceReal.values(), outfile)

            binLabel = 'qcd' # draw inputs, num and den
            frameTitle = 'qcd elec: '+reg+' tight;'+var
            histoName = 'elec_'+reg+'_all_flavor_'+var+'_num'
            canvasName = 'elec_qcd'+reg+'_'+var+'_num'
            histosFlavorSliceQcd = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSliceQcd, canvasName, outplotdir, frameTitle)

            frameTitle = 'qcd elec: '+reg+' loose;'+var
            histoName = 'elec_'+reg+'_all_flavor_'+var+'_den'
            canvasName = 'elec_qcd'+reg+'_'+var+'_den'
            histosFlavorSliceQcd = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSliceQcd, canvasName, outplotdir, frameTitle)

            binLabel = 'conv' # draw inputs, num and den
            frameTitle = 'conv elec: '+reg+' tight;'+var
            histoName = 'elec_'+reg+'_all_flavor_'+var+'_num'
            canvasName = 'elec_conv'+reg+'_'+var+'_num'
            histosFlavorSliceConv = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSliceConv, canvasName, outplotdir, frameTitle)

            frameTitle = 'conv elec: '+reg+' loose;'+var
            histoName = 'elec_'+reg+'_all_flavor_'+var+'_den'
            canvasName = 'elec_conv'+reg+'_'+var+'_den'
            histosFlavorSliceConv = sliceHistosLeptonType(inputFiles, histoName, binLabel)
            plotStackedHistos(histosFlavorSliceConv, canvasName, outplotdir, frameTitle)

            # normalize and draw fractions (den only)
            histosFlavorSlice = dict([(k+'_qcd',  h) for k,h in histosFlavorSliceQcd.iteritems()] +
                                     [(k+'_conv', h) for k,h in histosFlavorSliceConv.iteritems()])
            normalizeHistos(histosFlavorSlice)
            frameTitle = 'qcd elec: '+reg+';'+var
            canvasBaseName = 'elec_qcd'+reg+'_'+var+'_frac'
            plotFractionsHistos (histosFlavorSliceQcd, canvasBaseName+'_histo', outplotdir, frameTitle)
            plotFractionsStacked(histosFlavorSliceQcd, canvasBaseName+'_stack', outplotdir, frameTitle)
            frameTitle = 'conv elec: '+reg+';'+var
            canvasBaseName = 'elec_conv'+reg+'_'+var+'_frac'
            plotFractionsHistos (histosFlavorSliceConv, canvasBaseName+'_histo', outplotdir, frameTitle)
            plotFractionsStacked(histosFlavorSliceConv, canvasBaseName+'_stack', outplotdir, frameTitle)
            frameTitle = 'elec: '+reg+';'+var
            canvasBaseName = 'elec_fake'+reg+'_'+var+'_frac'            
            plotFractionsStacked(histosFlavorSliceConv, canvasBaseName+'_stack', outplotdir, frameTitle,
                                 histos2=histosFlavorSliceQcd, h1Label='conv', h2Label='qcd')
            saveHistos(histosFlavorSliceQcd.values(), outfile)
            saveHistos(histosFlavorSliceConv.values(), outfile)



def sliceHistosLeptonType(inputFiles, histoName, binLabel) :
    "Histos have one lepton source for each x bin; take a slice for one specific bin"
    histos = dict((p, f.Get(histoName)) for p, f in inputFiles.iteritems())
    assert all(h for h in histos.values()),"missing histo '%s'\n%s"%(histoName, '\n'.join("%s : %s"%(k, str(v)) for k, v in histos.iteritems()))
    bin = first(histos).GetXaxis().FindBin(binLabel)
    actualBinLabels = dict((h.GetName(), h.GetXaxis().GetBinLabel(bin)) for h in histos.values())
    assert all(binLabel == v for h, v in actualBinLabels.iteritems()),"expecting bin '%s', got \n%s"%(binLabel, actualBinLabels)
    histosFlavorSlice = dict((p, h.ProjectionY(p+'_'+h.GetName()+'_'+binLabel, bin, bin, 'e')) for p, h in histos.iteritems())
    return histosFlavorSlice
def guessBaseHistoname(histonames=[]) :
    basename = commonSuffix(histonames)
    basename = commonPrefix(histonames) if not len(basename) else basename
    basename = histonames[0] if not len(basename) else basename
    return basename.strip(' _')
def normalizeHistos(histosFlavorSlice) :
    basename = guessBaseHistoname([h.GetName() for h in histosFlavorSlice.values()])
    tot = first(histosFlavorSlice).Clone(basename+'_tot')
    tot.Reset()
    for h in histosFlavorSlice.values() : tot.Add(h)
    for b in getBinIndices(tot) : tot.SetBinError(b, 0.0) # norm is a constraint, without error
    for h in histosFlavorSlice.values() : h.Divide(tot)
def saveHistos(histos=[], outfile=None) :
    outfile.cd()
    for h in histos : h.Write()
def plotStackedHistos(histosFlavorSlice={}, canvasName='', outputDir='./', frameTitle='stack') :
    "Plot the input histos used to compute the fractions"
    histos = histosFlavorSlice
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    stack = r.THStack('stack_'+canvasName, '')
    leg = topRightLegend(can, 0.275, 0.475, shift=-0.025)
    leg.SetBorderSize(0)
    colors = SampleUtils.colors
    procs = sorted(histos.keys())
    for p in procs:
        h = histos[p]
        h.SetFillColor(colors[p])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
        stack.Add(h)
    for s in procs[::-1] : leg.AddEntry(histos[s], s, 'F') # stack goes b-t, legend goes t-b
    stack.Draw('hist e')
    leg.Draw('same')
    tex = r.TLatex()
    tex.SetNDC(True)
    tex.SetTextFont(first(histos).GetTitleFont())
    tex.SetTextSize(first(histos).GetTitleSize())
    tex.DrawLatex(0.1, 0.925, frameTitle.split(';')[0])
    can.Update() # force stack to create canMaster
    canMaster = stack.GetHistogram()
    canMaster.SetTitle(frameTitle)
    canMaster.Draw('axis same')
    can._graphical_objects = [stack, canMaster, leg, tex] + [h for h in stack.GetStack()]
    can.Update()
    for ext in ['png','eps'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)

def plotFractionsHistos(histos={}, canvasName='', outputDir='./', frameTitle='title;p_{T} [GeV]; efficiency', zoomIn=False) :
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    padMaster = None
    colors, markers = SampleUtils.colors, SampleUtils.markers
    for s,h in histos.iteritems() :
        h.SetLineColor(colors[s] if s in colors else r.kBlack)
        h.SetMarkerColor(h.GetLineColor())
        h.SetMarkerStyle(markers[s] if s in markers else r.kFullCircle)
        drawOpt = 'ep same' if padMaster else 'ep'
        h.Draw(drawOpt)
        if not padMaster : padMaster = h
    minY, maxY = getMinMax(histos.values()) if zoomIn else (0.0, 1.0)
    padMaster.GetYaxis().SetRangeUser(min([0.0, minY]), 1.1*maxY)
    padMaster.SetTitle(frameTitle)
    padMaster.SetStats(False)
    drawLegendWithDictKeys(can, histos)
    can.Update()
    for ext in ['png','eps'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)

def plotFractionsStacked(histos={}, canvasName='', outputDir='./', frameTitle='title;p_{T} [GeV]; efficiency',
                         histos2={}, h1Label='', h2Label='') :
    "plot stack of histos; optional histos2 are also stacked with the same colors, but with different fill style"
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    stack = r.THStack('stack_'+canvasName, '')
    leg = topRightLegend(can, 0.275, 0.475, shift=+0.010)
    leg.SetBorderSize(0)
    colors = SampleUtils.colors
    procs = sorted(histos.keys())
    def setHistoAtts(histo, process, colors=colors, fillStyle=None) :
        h, p = histo, process
        h.SetMarkerSize(0)
        h.SetFillColor(colors[p])
        h.SetLineColor(h.GetFillColor())
        h.SetDrawOption('bar')
        if fillStyle : h.SetFillStyle(fillStyle)
    for p in procs:
        h = histos[p]
        setHistoAtts(h, p)
        stack.Add(h)
    for s in procs[::-1] : leg.AddEntry(histos[s], s, 'F') # stack goes b-t, legend goes t-b
    for p in filter(lambda p: p in histos2, procs) :
        setHistoAtts(histos2[p], p, fillStyle=3001)
        stack.Add(histos2[p])
    stack.Draw('hist')
    if h1Label or h2Label:
        leg.SetHeader(', '.join(["%s : %s"%(f,l) if l else '' for f,l in [('darker',h1Label), ('lighter',h2Label)]]))
    leg.Draw('same')
    tex = r.TLatex()
    tex.SetNDC(True)
    tex.DrawLatex(0.1, 0.9, 'fractions: '+frameTitle[:frameTitle.find(';')-1])
    can.Update() # force stack to create canMaster
    canMaster = stack.GetHistogram()
    canMaster.SetTitle(frameTitle)
    canMaster.Draw('axis same')
    can._graphical_objects = [stack, canMaster, leg, tex] + [h for h in stack.GetStack()]
    can.Update()
    for ext in ['png','eps'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)

if __name__=='__main__' :
    main()
