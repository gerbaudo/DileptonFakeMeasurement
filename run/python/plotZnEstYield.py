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

import collections, optparse, sys
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

parser = optparse.OptionParser()
parser.add_option("-S", "--scale-sig", dest="sigScale", default=defaultSigScale, type='float',
                  help="scale the signal yield by this factor (default %.1f)" % defaultSigScale)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                  help="print more details about what is going on")
(options, args) = parser.parse_args()
sigScale        = options.sigScale
verbose         = options.verbose

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
    templateLine = 'Dataset MC1,MN2[GeV]  MN1[GeV]' # + m selections
    minNfields = len(templateLine.split())
    lines = open(filename).readlines()
    fields = lines[0].split()
    assert len(fields) > minNfields,"not enough columns (%d), should be >%d"%(len(fields), minNfields)
    dsIdx, mc1Idx, mn1Idx = fields.index('Dataset'), fields.index('MC1,MN2[GeV]'), fields.index('MN1[GeV]')
    assert dsIdx, 'cannot find ds index'
    assert mc1Idx, 'cannot find mc1 index'
    assert mn1Idx, 'cannot find mn1 index'
    selections = line0.split()[3:] # the first 3 columns are not yields
    allCounts = {}
    for line in lines[1:] :
        line = line.strip()
        fields = line.split()
        ds, mc1, mn1 = fields[dsIdx], fields[mc1Idx], fields[mn1Idx]
        counts = [float(c) for c in fields[3:]]
        key = {'ds':ds, 'mc1':mc1, 'mn1':mn1}
        assert key not in allCounts, "duplicate entry %s"%str(key)
        allCounts[key] = dict(zip(selections, counts))
    return allCounts
