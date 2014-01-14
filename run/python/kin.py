
# kinematic utilities
#
# davide.gerbaudo@gmail.com
# Jan 2014

import math

import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options

fabs, cos, sin, pi, sqrt = math.fabs, math.cos, math.sin, math.pi, math.sqrt
def phi_mpi_pi(phi) :
    pi = math.pi
    while phi < -pi : phi += pi
    while phi > +pi : phi -= pi
    return phi

tlv = r.TLorentzVector
def FourMom2TLorentzVector(fm) :
    l = tlv()
    l.SetPxPyPzE(fm.px, fm.py, fm.pz, fm.E)
    return l
def addTlv(l) :
    l.p4 = FourMom2TLorentzVector(l)
    return l
def computeMt(lep, met) :
    return sqrt(2.0 * lep.Pt() * met.Et() *(1.0-cos(lep.DeltaPhi(met))))
def computeHt(met, leptsJets=[]) :
    return sum(o.Pt() for o in leptsJets) + met.Et()
def computeMetRel(met, leptsJets=[]) :
    minDphi = min([0.5*pi]+[fabs(met.DeltaPhi(o)) for o in leptsJets])
    return met.Et()*sin(minDphi)
def lepIsSeparatedFromOther(l, otherLeps=[], minDr=0.05) :
    return all(l.DeltaR(ol)>minDr for ol in otherLeps)
def lepPairIsZcand(l0, l1) :
    "Inputs are susy::wh::FourMom"
    def lepFlavor(l) : return 'el' if l.isEl else 'mu' if l.isMu else 'other'
    l0Fl, l1Fl = lepFlavor(l0), lepFlavor(l1)
    elOrMu = l0Fl in ['el','mu']
    sameFlavor = l0Fl==l1Fl
    oppCharge  = l0.charge*l1.charge < 0.0
    return elOrMu and sameFlavor and oppCharge
def deltaMZ0((la, lb)) :
    "given a pair of leptons, return the abs difference m_ll - m_Z"
    return abs((la + lb).M() - 91.2)
def getDilepType(fmLep0, fmLep1) :
    "Given two susy::wh::FourMom, return ee/em/mm"
    def FourMom2LepType(fm) : return 'e' if fm.isEl else 'm' if fm.isMu else None
    dilepType = ''.join(sorted(FourMom2LepType(l) for l in [fmLep0, fmLep1])) # sort -> em instead of me
    return dilepType
