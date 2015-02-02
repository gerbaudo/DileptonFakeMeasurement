
import numpy as np

from utils import first, rmIfExists
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
def mdeltarBinEdges() : return np.array([0.0, 20.0, 40.0, 60.0, 100.0, 200.0])

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
    padMaster.SetMaximum(0.25 if (maxY < 0.25 and zoomIn) else 1.1)
    padMaster.SetTitle(frameTitle)
    padMaster.SetStats(False)
    drawLegendWithDictKeys(can, effs)
    can.Update()
    for ext in ['png','eps'] :
        outFilename = outputDir+'/'+canvasName+'.'+ext
        rmIfExists(outFilename)
        can.SaveAs(outFilename)

def plot2dEfficiencies(effs={}, canvasName='', outputDir='./', frameTitle='efficiency; #eta; p_{T} [GeV]', textFormat='.3f', zoomIn=False) :
    can = r.TCanvas(canvasName, '', 800, 600)
    can.cd()
    origTextFormat = r.gStyle.GetPaintTextFormat()
    r.gStyle.SetPaintTextFormat(textFormat)
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

#___________________________________________________________
def isTight_std(l):
    "tight with standard isolation"
    return elecIsTight_std  (l) if l.isEl else muonIsTight_std  (l) if l.isMu else False
def isTight_tight(l):
    "tight with standard isolation, but tighter cuts"
    return elecIsTight_tight(l) if l.isEl else muonIsTight_tight(l) if l.isMu else False
def isTight_minden(l):
    "tight with modified isolation: denominator=min(pt,60) rather than pt"
    return elecIsTight_minden(l) if l.isEl else muonIsTight_minden(l) if l.isMu else False
def isTight_wh (l):
    "tight with tighter cuts and modified isolation"
    return elecIsTight_wh   (l) if l.isEl else muonIsTight_wh   (l) if l.isMu else False

def elecIsTight_std   (l): return l.isTightPp and elecIsFromPv(l) and elecIsIsolated(l, denominator_std(l), etConeThres=0.18, ptConeThres=0.16)
def elecIsTight_tight (l): return l.isTightPp and elecIsFromPv(l) and elecIsIsolated(l, denominator_std(l), etConeThres=0.13, ptConeThres=0.07)
def elecIsTight_minden(l): return l.isTightPp and elecIsFromPv(l) and elecIsIsolated(l, denominator_wh (l), etConeThres=0.18, ptConeThres=0.16)
def elecIsTight_wh    (l): return l.isTightPp and elecIsFromPv(l) and elecIsIsolated(l, denominator_wh (l), etConeThres=0.13, ptConeThres=0.07)

def muonIsTight_std   (l): return muonIsFromPv(l) and muonIsIsolated(l, denominator_std(l), etConeThres=None, ptConeThres=0.12)
def muonIsTight_tight (l): return muonIsFromPv(l) and muonIsIsolated(l, denominator_std(l), etConeThres=0.14, ptConeThres=0.06)
def muonIsTight_minden(l): return muonIsFromPv(l) and muonIsIsolated(l, denominator_wh(l),  etConeThres=None, ptConeThres=0.12)
def muonIsTight_wh    (l): return muonIsFromPv(l) and muonIsIsolated(l, denominator_wh (l), etConeThres=0.14, ptConeThres=0.06)

def lepIsTight_std   (l): return muonIsTight_std   (l) if l.isMu else elecIsTight_std   (l)
def lepIsTight_tight (l): return muonIsTight_tight (l) if l.isMu else elecIsTight_tight (l)
def lepIsTight_minden(l): return muonIsTight_minden(l) if l.isMu else elecIsTight_minden(l)
def lepIsTight_wh    (l): return muonIsTight_wh    (l) if l.isMu else elecIsTight_wh    (l)


# isolation from Liz (see email 2014-06-23, subj 'Fake rate for SS EWK')
def lepIsTight_07(l):
    return (muonIsFromPv(l) and
            muonIsIsolated(l, denominator_wh(l), etConeThres=0.14, ptConeThres=0.06)
            if l.isMu else
            l.isTightPp and
            elecIsFromPv(l) and
            elecIsIsolated(l, denominator_wh(l), etConeThres=0.13, ptConeThres=0.07))
def lepIsTight_06(l):
    return (muonIsFromPv(l) and
            muonIsIsolated(l, denominator_wh(l), etConeThres=0.12, ptConeThres=0.06)
            if l.isMu else
            l.isTightPp and
            elecIsFromPv(l) and
            elecIsIsolated(l, denominator_wh(l), etConeThres=0.12, ptConeThres=0.06))
def lepIsTight_05(l):
    return (muonIsFromPv(l) and
            muonIsIsolated(l, denominator_wh(l), etConeThres=0.11, ptConeThres=0.05)
            if l.isMu else
            l.isTightPp and
            elecIsFromPv(l) and
            elecIsIsolated(l, denominator_wh(l), etConeThres=0.11, ptConeThres=0.05))



def denominator_wh (l) :
    pt = l.p4.Pt()
    return 1.0/min([pt, 60.0]) if pt>0.0 else None

def denominator_std(l) :
    pt = l.p4.Pt()
    return 1.0/pt if pt>0.0 else None

def elecIsFromPv(l):
    # maxD0Sig, maxZ0SinTheta = 3.0, 0.4 # see SusyNtTools::isSignalElectron, WH values
    maxD0Sig, maxZ0SinTheta = 5.0, 0.4 # see SusyNtTools::isSignalElectron, std values
    return abs(l.d0Signif) < maxD0Sig and abs(l.z0SinTheta) < maxZ0SinTheta

def muonIsFromPv(l):
    maxD0Sig, maxZ0SinTheta = 3.0, 1.0 # see SusyNtTools::isSignalMuon
    return abs(l.d0Signif) < maxD0Sig and abs(l.z0SinTheta) < maxZ0SinTheta

def lepIsFromPv(l):
    return elecIsFromPv(l) if l.isEl else muonIsFromPv(l) if l.isMu else False

def elecIsIsolated(l, denom, etConeThres, ptConeThres):
    pt, etCone, ptCone = l.p4.Pt(), l.etConeCorr, l.ptConeCorr
    return ((etCone*denom < etConeThres if etConeThres else True) and
            (ptCone*denom < ptConeThres if ptConeThres else True)
            if denom else False)

def muonIsIsolated(l, denom, etConeThres, ptConeThres):
    pt, etCone, ptCone = l.p4.Pt(), l.etConeCorr, l.ptConeCorr # DG: Nt uses ptConeCorr only for 3lep?
    return ((etCone*denom < etConeThres if etConeThres else True) and
            (ptCone*denom < ptConeThres if ptConeThres else True)
            if denom else False)

def fetchSfHistos(inputSfFiles=[], lepton='', verbose=False):
    from compute_fake_el_scale_factor import histoname_sf_vs_eta
    fileNames = inputSfFiles
    histos = dict()
    print "fetchSfHistos: fileNames ",fileNames
    if not (type(fileNames)==list and
            len(fileNames) in [1, 2]):
        print "fetchSfHistos expects one or two files (hflf+conv), got %s"%str(inputSfFiles)
        print "returning ",histos
        return histos
    if verbose : print "retrieving scale factors from %s"%inputSfFiles
    fname_hflf = first(filter(lambda _ : 'hflf' in _, fileNames))
    fname_conv = first(filter(lambda _ : 'conv' in _, fileNames))
    file_hflf = r.TFile.Open(fname_hflf) if fname_hflf else None
    file_conv = r.TFile.Open(fname_conv) if fname_conv else None
    hname = histoname_sf_vs_eta(lepton)
    histo_hflf = file_hflf.Get(hname) if file_hflf else None
    histo_conv = file_conv.Get(hname) if file_conv else None
    if histo_hflf : histos['hflf'] = composeEtaHistosAs2dPtEta(input1Dhisto=histo_hflf,
                                                               outhistoname=hname+'_hflf')
    if histo_conv : histos['conv'] = composeEtaHistosAs2dPtEta(input1Dhisto=histo_conv,
                                                               outhistoname=hname+'_conv')
    for f in [file_hflf, file_conv] :
        if f : f.Close()
    return histos

def composeEtaHistosAs2dPtEta(input1Dhisto=None, outhistoname='') :
    """take the 1D scale factor histogram (vs eta), and build a 2D
    histo that has (pt,eta) on (x,y); see
    MeasureFakeRate2::initHistos"""
    ptBins = ptBinEdges()
    etaBins = etaBinEdges()
    h = r.TH2F(outhistoname, '', len(ptBins)-1, ptBins, len(etaBins)-1, etaBins)
    h.SetDirectory(0)
    h.Sumw2()
    for iX in range(1, 1+len(ptBins)):
        for iY in range(1, 1+len(etaBins)):
            h.SetBinContent(iX, iY, input1Dhisto.GetBinContent(iY))
            h.SetBinError  (iX, iY, input1Dhisto.GetBinError  (iY))
    return h

#___________________________________________________________
# see SusyDefs.h:TrigBit
triggerBitNames = ['e7_medium1', #2012 triggers
                   'e12Tvh_loose1',
                   'e12Tvh_medium1',
                   'e24vh_medium1',
                   'e24vhi_medium1',
                   '2e12Tvh_loose1',
                   'e24vh_medium1_e7_medium1',
                   'mu8',
                   'mu13',
                   'mu18_tight',
                   'mu24i_tight',
                   '2mu13',
                   'mu18_tight_mu8_EFFS',
                   'e12Tvh_medium1_mu8',
                   'mu18_tight_e7_medium1',
                   'g20_loose', # Photon Triggers
                   'g40_loose',
                   'g60_loose',
                   'g80_loose',
                   'g100_loose',
                   'g120_loose',
                   'tau20_medium1', # Tau triggers
                   'tau20Ti_medium1',
                   'tau29Ti_medium1',
                   'tau29Ti_medium1_tau20Ti_medium1',
                   'tau20Ti_medium1_e18vh_medium1',
                   'tau20_medium1_mu15',
                   'e18vh_medium1', # Missing trigger flags for lep-tau matching
                   'mu15',
                   '2mu8_EFxe40wMu_tclcw', # MissingEt trigger
                   ]
def triggerBit(bitName): return 1<<triggerBitNames.index(bitName)
def passTrigger(l, trigname): return l.trigFlags & triggerBit(trigname)
#___________________________________________________________
tupleStemsAndNames = [ # stem filename, treename (see MeasureFakeRate2::Begin)
    ("hflf"     , "HeavyFlavorControlRegion"),
    ("hflfss"   , "HeavyFlavorSsControlRegion"),
    ("conv"     , "ConversionControlRegion"),
    ("zmmejets" , "ZmmeVetoPlusJetsRegion"),
    ("ssinc"    , "SameSignRegion"),
    ("ssinc1j"  , "SameSign1jetControlRegion"),
    ("mcconv"   , "ConversionExtractionRegion"),
    ("mcqcd"    , "HfLfExtractionRegion"),
    ("mcreal"   , "RealExtractionRegion"),
    ("emu"      , "EmuRegion"),
    ("razor0j"  , "Razor0jRegion"),
    ("razor1j"  , "Razor1jRegion"),
    ]
#___________________________________________________________
