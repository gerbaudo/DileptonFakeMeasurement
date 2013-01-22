#!/bin/env python

# Utility classes and functions to manage samples and subsamples
#
# davide.gerbaudo@gmail.com
# Jan 2013

import ROOT as r

colors = {
    'ttbar'     : r.kRed+1,
    'zjets'     : r.kOrange-2,
    'wjets'     : r.kBlue-2,
    'diboson'   : r.kSpring+1,
    'singletop' : r.kAzure-4,
    'multijet'  : r.kWhite
    }

def guessSampleFromFilename(filename='', verbose=False) :
    if 'top_' in filename : return 'ttbar'
    elif 'Zjet_' in filename : return 'zjets'
    elif 'ZZ_' in filename \
         or 'WW_' in filename \
         or 'WZ_' in filename : return 'diboson'
    elif 'Wjet_' in filename : return 'wjets'
    else :
        if verbose : print "cannot guess samplename for %s" % filename
