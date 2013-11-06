#!/bin/env python
import sys
import ROOT as r
r.gROOT.SetBatch(True)

input  = r.TFile(sys.argv[1])
knames = sorted([k.GetName() for k in input.GetListOfKeys()])
print '\n'.join(knames)
