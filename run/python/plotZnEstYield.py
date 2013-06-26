#!/bin/env python

# Starting from the yield estimates from Anyes, produce Zn plots
#
# Inputs:
# - txt files(for sig and bkg) with estimated yields from grouped efficiencies
# Steps:
# - parse files for ee/em/mm
# - for each signal sample, compute Zn for each selection and pick the
#   selection providing the best Zn
# - combine ee/em/mm by doing a sum in quadrature of the Zn values (PWC)
# davide.gerbaudo@gmail.com
# Jun 2013

import collections
import glob
import math
import operator
import optparse
import sys
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from PickleUtils import readFromPickle
from SampleUtils import ModeAWhDbPar, ModeAWhDbReqid

#########
# default parameters [begin]
defaultSigScale    = 1.0
# default parameters [end]
#########

def parseBkgYields(filename) :
    "parse a file with 2 lines: header and counts"
    lines = open(filename).readlines()
    assert len(lines)==2, "cannot parse %s, %d lines"%(filename, len(lines))
    selections = lines[0].strip().split()
    counts = lines[1].split()[1:] # the first column is the label 'BKG'
    assert len(counts)==len(selections),"got %d counts and %s selections"%(len(counts), len(selections))
    counts = [float(c) for c in counts]
    return dict(zip(selections, counts))
def parseSigYields(filename) :
    """Parse a file with n lines + 1 header line.
    Each line has 3 parameters and m(>=1) counts for m selections.
    """
    templateLine = 'Dataset "MC1,MN2[GeV]"  MN1[GeV]' # + m selections; note xls puts some "
    minNfields = len(templateLine.split())
    lines = open(filename).readlines()
    fields = lines[0].split()
    assert len(fields) > minNfields,"not enough columns (%d), should be >%d"%(len(fields), minNfields)
    dsIdx, mc1Idx, mn1Idx = fields.index('Dataset'), fields.index('"MC1,MN2[GeV]"'), fields.index('MN1[GeV]')
    selections = fields[3:] # the first 3 columns are not yields
    allCounts = {}
    for line in lines[1:] :
        line = line.strip()
        if not line : continue # skip empty lines
        fields = line.split()
        ds, mc1, mn1 = fields[dsIdx], fields[mc1Idx], fields[mn1Idx]
        counts = [float(c) for c in fields[3:]]
        key = frozenset({'ds':ds, 'mc1':mc1, 'mn1':mn1}.items()) # must be hashable
        assert key not in allCounts, "duplicate entry %s"%str(key)
        allCounts[key] = dict(zip(selections, counts))
    return allCounts
def channelFromFilename(filename) :
    filename = filename.lower()
    if   '_ee' in filename : return 'ee'
    elif '_em' in filename : return 'em'
    elif '_mm' in filename : return 'mm'
def sigBkgFromFilename(filename) :
    filename = filename.lower()
    if   '_signal' in filename : return 'sig'
    elif '_bkg'    in filename : return 'bkg'
def findBestZn(countsPerSelBkg={}, countsPerSelSig={},
               sigScale=1.0, bkgRelErr=0.2) :
    selections = countsPerSelSig.keys()
    znFunc = r.RooStats.NumberCountingUtils.BinomialExpZ
    zns = dict([(sel, znFunc(sigScale*countsPerSelSig[sel],
                             countsPerSelBkg[sel],
                             bkgRelErr))
                for sel in selections])
    return max(zns.iteritems(), key=operator.itemgetter(1)) # returns (sel, Zn_value)
def buildPadMaster(points, histoname='padmaster', histotitle='') :
    points = [dict(p) for p in points]
    xLabel, yLabel = 'mc1', 'mn1'
    xs, ys = [p[xLabel] for p in points], [p[yLabel] for p in points]
    xRange = {'min': min(xs), 'max' : max(xs)}
    yRange = {'min': min(ys), 'max' : max(ys)}
    return r.TH2F(histoname,
                  histotitle+';'+xLabel+';'+yLabel,
                  50, float(xRange['min']), float(xRange['max']),
                  50, float(yRange['min']), float(yRange['max']))
def plotBestSelection(pointsWithSel, histoname, histotitle) :
    pm = buildPadMaster(pointsWithSel, histoname, histotitle)
    points = [(p['mc1'], p['mn1'], p['sel']) for p in [dict(pp) for pp in pointsWithSel]]
    xs, ys, sels = zip(*points) # vertical slice
    colors = [r.kRed, r.kBlue, r.kOrange, r.kViolet]
    uniqMarkers = dict()
    for i,s in enumerate(set(sels)) :
        m = r.TMarker(0.0, 0.0, r.kFullCircle)
        m.SetMarkerColor(colors[i])
        uniqMarkers[s] = m
    c = r.TCanvas('c_'+histoname, histoname, 800, 600)
    c.cd()
    pm.SetStats(0)
    pm.Draw('axis')
    allMarkers = []
    for x,y,sel in points :
        m = r.TMarker(float(x), float(y), r.kFullCircle)
        m.SetMarkerColor(uniqMarkers[sel].GetMarkerColor())
        m.Draw()
        allMarkers.append(m)
    leg = r.TLegend(0.1, 0.75, 0.4, 0.9, 'Opt. selection')
    leg.SetFillColor(0)
    leg.SetBorderSize(1)
    leg.Draw()
    for s,m in uniqMarkers.iteritems() : leg.AddEntry(m, s, 'p')
    c.Update()
    for ext in extensions : c.SaveAs(c.GetName()+'.'+ext)

def plotBestZns(znDict, histoname, histotitle, alsoNegative=False) :
    pm = buildPadMaster(znDict.keys(), histoname, histotitle)
    points = [(p['mc1'], p['mn1'], zn) for p, zn in [(dict(pp), v) for pp,v in znDict.iteritems()]]
    gr = r.TGraph2D(1)
    gr.SetTitle('g_'+histoname)
    gr.SetMarkerStyle(r.kFullSquare)
    gr.SetMarkerSize(2*gr.GetMarkerSize())
    for i, (x, y, zn) in enumerate(points) :
        if zn>0.0 or alsoNegative : gr.SetPoint(i+1, float(x), float(y), zn)
    c = r.TCanvas('c_'+histoname, histoname, 800, 600)
    c.cd()
    pm.SetStats(0)
    pm.Draw('axis')
    gr.Draw('colz same')
    tex = r.TLatex(0.0, 0.0, '')
    tex.SetTextFont(pm.GetTitleFont())
    for x, y, zn in points : tex.DrawLatex(float(x), float(y), "%.2f"%zn)
    c.Update()
    for ext in extensions : c.SaveAs(c.GetName()+'.'+ext)

def combineZn(znDictPerChannel={}, alsoNegative=False) :
    """combine the zn from multiple channels by summing them in
    quadrature; return a dict with the same keys(except the selection)
    as input and val=combined zn
    """
    sqrt = math.sqrt
    values = collections.defaultdict(list)
    for ch, valsPerCh in znDictPerChannel.iteritems() :
        for k, v in valsPerCh.iteritems() :
            k = dict(k)
            k = frozenset(dict([(p, k[p]) for p in ['ds','mc1','mn1']]).items()) # don't care about sel
            values[k].append(v)
    return dict([(k, sqrt(sum([v*v for v in vv if v>0.0 or alsoNegative]))) for k, vv in values.iteritems()])

#_______________________________________

if __name__=='__main__' :
    parser = optparse.OptionParser()
    parser.add_option("-S", "--scale-sig", dest="sigScale", default=defaultSigScale, type='float',
                      help="scale the signal yield by this factor (default %.1f)" % defaultSigScale)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                      help="print more details about what is going on")
    (options, args) = parser.parse_args()
    sigScale        = options.sigScale
    verbose         = options.verbose
    extensions      = ['eps','png']
    setupLabel = "(SS, signal x%.1f)"%sigScale

    filesBkg = glob.glob('estimatedYield/*bkg*txt')
    filesSig = glob.glob('estimatedYield/*signal*txt')

    filesPerChannel = collections.defaultdict(list)
    for f in filesSig + filesBkg : filesPerChannel[channelFromFilename(f)].append(f)
    def twoFilenamesAreValid(filenames) :
        'exactly one sig and one bkg'
        f0, f1 = filenames[0], filenames[1]
        type0, type1 = sigBkgFromFilename(f0), sigBkgFromFilename(f1)
        return len(filenames)==2 and type0 and type1 and type0 != type1
    assert all([twoFilenamesAreValid(ff) for ff in filesPerChannel.values()]),"all channels should have 1sig + 1bkg file"

    optimalZnPerChannel = dict()
    for ch, twoFiles in filesPerChannel.iteritems() :
        fnameSig, fnameBkg = twoFiles
        if sigBkgFromFilename(fnameSig) is not 'sig' : fnameSig, fnameBkg = fnameBkg, fnameSig
        countsBkg  = parseBkgYields(fnameBkg)
        countsSigs = parseSigYields(fnameSig)
        znThisChannel = dict()
        for sig, countsSig in countsSigs.iteritems() :
            bestSel, bestZn = findBestZn(countsBkg, countsSig, sigScale)
            key = dict(sig)
            key['sel'] = bestSel # add a keyval pair
            key = frozenset(key.items())
            assert key not in znThisChannel,"duplicate optimal selection? %s"%str(key)
            znThisChannel[key] = bestZn
        plotBestZns(znThisChannel, ch+'_optZn', 'Best Zn: '+ch+' '+setupLabel)
        plotBestSelection(znThisChannel.keys(), ch+'_optSel', 'Optimal selection: '+ch+' '+setupLabel)
        optimalZnPerChannel[ch] = znThisChannel
    optimalZnPerChannel['ll'] = combineZn(optimalZnPerChannel)
    plotBestZns(optimalZnPerChannel['ll'], 'll_optZn', 'Best Zn: ll '+setupLabel)
