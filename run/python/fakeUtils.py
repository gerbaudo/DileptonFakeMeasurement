
import numpy as np

from utils import rmIfExists
from rootUtils import importRoot, buildRatioHistogram, getMinMax, drawLegendWithDictKeys
r = importRoot()
import SampleUtils

def samples() : return ['allBkg', 'ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def fakeProcesses() : return ['ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def getInputFiles(inputDirname, tag, verbose=False) :
    inDir = inputDirname
    tag = tag if tag.startswith('_') else '_'+tag
    files = dict(zip(samples(), [r.TFile.Open(inDir+'/'+s+tag+'.root') for s in samples()]))
    if verbose : print "getInputFiles('%s'):\n\t%s"%(inputDirname, '\n\t'.join("%s : %s"%(k, f.GetName()) for k, f in files.iteritems()))
    return files
def buildRatio(inputFile=None, histoBaseName='') :
    num, den = inputFile.Get(histoBaseName+'_num'), inputFile.Get(histoBaseName+'_den')
    return buildRatioHistogram(num, den, histoBaseName +'_rat')

def mtBinEdges() : return np.array([0.0, 20.0, 40.0, 60.0, 100.0, 200.0])
def ptBinEdges() : return np.array([10.0, 20.0, 35.0, 100.0])
def etaBinEdges() : return np.array([0.0, 1.37, 2.50])

def leptonTypes() : return ['tight', 'loose', 'real_tight', 'real_loose', 'fake_tight', 'fake_loose']
def allLeptonSources() : return ['heavy',   'light', 'conv',  'real',  'qcd', 'unknown'] # see FakeLeptonSources.h
def enum2source(l):
    "convert the int(enum) stored in the tree to a 'lepton source' string"
    return allLeptonSources()[l.source]
def leptonSources() : return [s for s in allLeptonSources() if s not in ['qcd']] # qcd is just hf+lf
def colorsFillSources() :
    return dict(zip(leptonSources(), [r.kBlue-10, r.kMagenta-10, r.kRed-8, r.kGreen-6, r.kCyan-6, r.kGray+1]))
def colorsLineSources() :
    return dict(zip(leptonSources(), [r.kBlue, r.kMagenta, r.kRed, r.kGreen, r.kCyan, r.kGray+1]))
def markersSources() :
    return dict(zip(leptonSources(), [r.kPlus, r.kCircle, r.kMultiply, r.kOpenSquare, r.kOpenTriangleUp, r.kOpenCross]))

def plot1dEfficiencies(effs={}, canvasName='', outputDir='./', frameTitle='title;p_{T} [GeV]; efficiency', zoomIn=False) :
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    padMaster = None
    colors, markers = SampleUtils.colors, SampleUtils.markers
    for s,h in effs.iteritems() :
        h.SetLineColor(colors[s] if s in colors else r.kBlack)
        h.SetMarkerColor(h.GetLineColor())
        h.SetMarkerStyle(markers[s] if s in markers else r.kFullCircle)
        drawOpt = 'ep same' if padMaster else 'ep'
        h.Draw(drawOpt)
        if not padMaster : padMaster = h
    minY, maxY = getMinMax(effs.values()) if zoomIn else (0.0, 1.0)
    padMaster.GetYaxis().SetRangeUser(min([0.0, minY]), 1.1*maxY)
    padMaster.SetMinimum(0.0)
    padMaster.SetMaximum(1.1*maxY)
    padMaster.SetMaximum(0.25)
    padMaster.SetTitle(frameTitle)
    padMaster.SetStats(False)
    drawLegendWithDictKeys(can, effs)
    can.Update()
    for ext in ['png','eps'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)

def plot2dEfficiencies(effs={}, canvasName='', outputDir='./', frameTitle='efficiency; #eta; p_{T} [GeV]', zoomIn=False) :
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    origTextFormat = r.gStyle.GetPaintTextFormat()
    r.gStyle.SetPaintTextFormat('.2f')
    for s,h in effs.iteritems() :
        can.Clear()
        # todo minZ, maxZ = getMinMax(effs.values()) if zoomIn else (0.0, 1.0)
        minZ, maxZ = (0.0, 1.0)
        h.SetMarkerSize(1.5*h.GetMarkerSize())
        h.Draw('colz')
        h.Draw('text e same')
        h.GetZaxis().SetRangeUser(min([0.0, minZ]), maxZ)
        def dropCrPrefix(sr) : return sr.replace('CR_', '')
        h.SetTitle(dropCrPrefix(s)+' : '+frameTitle)
        h.SetStats(False)
        can.Update()
        for ext in ['png','eps'] :
            outFilename = outputDir+'/'+canvasName+'_'+s+'.'+ext
            rmIfExists(outFilename)
            can.SaveAs(outFilename)
    r.gStyle.SetPaintTextFormat(origTextFormat)
