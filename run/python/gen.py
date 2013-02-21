# largely based on
# https://github.com/betchart/susycaf/blob/master/calculables/gen.py
import collections, ROOT as r
#from supy import utils,analysisStep
#####################################
try:
    import pdgLookup
    pdgLookupExists = True
except ImportError:
    pdgLookupExists = False
#####################################


class PdgLookup:
    def __init__(self) :
        self.names = {
            -1:'/d', 1:'d',
             -2:'/u', 2:'u',
             -3:'/s', 3:'s',
             -4:'/c', 4:'c',
             -5:'/b', 5:'b',
             -6:'/t', 6:'t',
             -11:'e+', 11:'e-',
             -12:'ve', 12:'ve',
             -13:'mu+', 13:'mu-',
             -14:'vmu', 14:'vmu',
             -15:'tau+', 15:'tau-',
             -16:'vtau', 16:'vtau',
             21:'g',
             -22:'gamma', 22:'gamma',
             -23:'Z', 23:'Z',
             -24:'W-', 24:'W+',
             -25:'h', 25:'h',
             -511:'/B0',511:'B0',
             -513:'/B*0',513:'B*0',
             -521:'B-',521:'B+',
             -523:'B*-',523:'B*+',
             -531:'/B0s',531:'B0s',
             -533:'/B*0s',533:'B*0s',
        }
    def pdgid_to_name(self,id) : return self.names[id] if id in self.names else 'unknown'

#-#####################################
#-class ParticleCountFilter(analysisStep) :
#-    def __init__(self, reqDict) :
#-        self.reqDict = reqDict
#-    def select (self,eventVars) :
#-        for key,value in self.reqDict.iteritems() :
#-            if eventVars["GenParticleCategoryCounts"][key]!=value : return False
#-        return True
#####################################
#####################################

MeV2GeV = 0.001

class particlePrinter(object) :

    def __init__(self,minPt=-1.0,minStatus=-1):
        self.minPt=minPt
        self.minStatus=minStatus
        
    def uponAcceptance (self,tree) :
        pdgLookupExists = True
        pdgLookup = PdgLookup()
        
        mc_parent_index = tree.mc_parent_index
        parents=set([p for pp in mc_parent_index for p in pp])
        
        #print "parents: ",parents
        print "-----------------------------------------------------------------------------------"
        print " i  st   par         id            name        pt       eta    phi    mass"
        print "-----------------------------------------------------------------------------------"

        size     = tree.mc_n
        pts      = tree.mc_pt
        etas     = tree.mc_eta
        phis     = tree.mc_phi
        pdgs     = tree.mc_pdgId
        masses   = tree.mc_m
        statuses = tree.mc_status
    
        maxPrintSize=75
        for iGen in range(min([maxPrintSize,size])) :
            pt = pts[iGen]
            if pt < self.minPt : continue
            status = statuses[iGen]
            if status<self.minStatus : continue

            pars = [i for i in mc_parent_index[iGen]]
            pdgId= pdgs[iGen]
            outString=""
            outString+="%#2d"%iGen
            outString+=" %#3d"%status
            outString+= str(pars).rjust(6)
            outString+=" %#10d"%pdgId
            if pdgLookupExists : outString+=" "+pdgLookup.pdgid_to_name(pdgId).rjust(15)
            else :                 outString+="".rjust(16)
            outString+="  %#8.1f"%(pts[iGen]*MeV2GeV)
            outString+="  %#8.1f"%etas[iGen]
            outString+="  %#5.1f"%phis[iGen]
            outString+="  %#6.1f"%(masses[iGen]*MeV2GeV)
            if not (iGen in parents) : outString+="   non-mo"
            print outString
        print

def findInterestingHiggsWithChiAndPar(pdg, parents, children) :
    #iHiggs = next((i for i,p in enumerate(pdg) if p==higgsPdg), None)
    higgsPdg = 25
    iHiggs = [i for i,p in enumerate(pdg) if p==higgsPdg]
    higgsChildren = [[pdg[i] if i<len(pdg) else 0 for i in hhc]
                     for hhc in [children[i]
                                 for i in iHiggs] ]
    higgsParents = [[pdg[i] if i<len(pdg) else 0 for i in hhp]
                    for hhp in [parents[i]
                                for i in iHiggs] ]
    # there can be several intermediate higgs
    boringIhiggs = [ih for i, ih in enumerate(iHiggs)
                    if len(higgsChildren[i])<2 or higgsPdg in higgsChildren[i]]
    interestingIhiggs = [i for i in iHiggs if i not in boringIhiggs]
    higgsChildren = [hc for hc,ih in zip(higgsChildren, iHiggs) if ih in interestingIhiggs]
    higgsParents = [hp for hp,ih in zip(higgsParents, iHiggs) if ih in interestingIhiggs]
    return interestingIhiggs, higgsChildren, higgsParents
def guessHdecayLabel(higgsChildren) :
    ch = higgsChildren
    if   any(p in ch for p in [-24, +24])  : return 'WW'
    elif any(p in ch for p in [23])        : return 'ZZ'
    elif any(p in ch for p in [-15, +15])  : return 'tautau'
    elif any(p in ch for p in [-5, +5])    : return 'bbbar'
    elif any(p in ch for p in [-13, +13])  : return 'mumu'
    else                                   : return 'unknown'
