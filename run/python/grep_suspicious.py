#!/bin/env python

# extract from the logs the printouts of the suspicious events in the fake control region
#
# Usage:
# grep_suspicious.py log/susyplot/*period*${TAG}.log
#
# davide.gerbaudo@gmail.com
# 2014-04-01

import sys

from rootUtils import importRoot
r = importRoot()
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

def main():
    electrons, muons, jets, mets = [], [], [], []
    for filename in sys.argv[1:]:
        with open(filename) as f:
            el, mu, jet, met = extractText(f)
            electrons += el
            muons += mu
            jets += jet
            mets += met
    for obj_name in ['electrons', 'muons', 'jets']:
        obj = eval(obj_name)
        plot_eta_phi(obj_name, obj)

def extractText(file):
    """Extract a block that looks like this:
    fake suspicious event run 210302 event 73750293 em 2j
    Run 210302 Event 73750293 Stream Muons w 1
    Mu : q  1 pt  46.31 eta  1.98 phi  2.69 cb 1 type 0 origin 0
    El : q  1 pt  10.47 eta  0.14 phi  1.66 tight 1 type 0 origin 0
    Met : pt  35.48 phi 0.65
    Jet : pt  59.82 eta  0.72 phi -1.12 jvf 0.97 mv1  0.05
    Jet : pt  27.30 eta  0.67 phi  2.44 jvf 0.47 mv1  0.06
    Jet : pt  21.93 eta -1.83 phi  1.34 jvf 0.07 mv1  0.27
    """
    inBlock = False
    def startBlock(l) : return l.startswith('fake suspicious event')
    def stillInBlock(l) : return any(l.startswith(t) for t in ['Run ', 'Mu :', 'El :', 'Met :', 'Jet :'])
    def floats_from_tokens(l, positions=[]) :
        tokens = l.strip().split()
        return (float(tokens[i]) for i in positions)
    def electron_from_line(l) : return floats_from_tokens(l, [7, 9])
    def muon_from_line(l) :     return floats_from_tokens(l, [7, 9])
    def jet_from_line(l) :      return floats_from_tokens(l, [5, 7])
    buffer = []
    electrons, muons, jets, mets = [], [], [], []
    for l in file.readlines():
        if not inBlock and startBlock(l) :
            print ''.join(buffer)
            inBlock = True
            buffer = [l]
        elif inBlock and stillInBlock(l) :
            buffer.append(l)
            if l.startswith('El :') : electrons.append(electron_from_line(l))
            elif l.startswith('Mu :') : muons.append(muon_from_line(l))
            elif l.startswith('Jet :') : jets.append(jet_from_line(l))
        else:
            inBlock = False
    if buffer:
        print ''.join(buffer)
    return electrons, muons, jets, mets
        
def plot_eta_phi(obj_name, objects):
    can = r.TCanvas(obj_name, obj_name+': eta phi', 800, 600)
    can.cd()
    h = r.TH2F(obj_name+'_eta_phi', obj_name, 60, -3.0, +3.0, 60, -3.14, +3.14)
    for x,y in objects:
        print x,y
        h.Fill(x, y, 1.0)
    h.Draw('col')
    can.SaveAs(obj_name+'.png')
    
if __name__=='__main__':
    main()
