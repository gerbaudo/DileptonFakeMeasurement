#!/bin/env python
#
# Script to dump the integral and raw counts just from one file
#
# ./printHistoYield.py file.root
#
# davide.gerbaudo@gmail.com
# 2013-06-05

import string
import sys
import ROOT as r

histos = {}
cw = colwidth = 8
for s in ["sr%d"%i for i in [6,7,8,9]] :
    histos[s] = ["%s_%s_onebin_NOM"%(s,l) for l in ['ee','em','mm']]

inputFile = r.TFile.Open(sys.argv[1])
field='%'+str(cw)+'s'

header = ''.join([field%s for s in [string.ljust('sr',cw),'ee','em','mm']])
print header
for s, hs in sorted(histos.iteritems()) :
    line  = string.ljust(s, cw)
    line += ''.join([field%("%.1f(%d)"%(h.Integral(), h.GetEntries()) if h.Integral() else "--")
                     for h in [inputFile.Get(h) for h in hs]])
    print line
print
