
import numpy as np

from rootUtils import importRoot, buildRatioHistogram
r = importRoot()

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
