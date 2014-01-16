#!/bin/env python

import math
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

def plot_zn_after_opt():
    znValues = znValuesAnyes
    znValues = znValuesDavide
    ds2p = convertDsidDict2points
    for ll in ['ee', 'mm', 'em'] :
        if ll not in znValues : continue
        for nj in ['eq1j', 'ge2j'] :
            if nj not in znValues[ll] : continue
            llnj = ll+'_'+nj
            plotZns(ds2p(znValues[ll][nj]),
                    histoname=llnj, histotitle=llnj+' ;mc_{1}; mn_{1}')
        plotZns(ds2p(combineZn(znValues[ll])),
                histoname=ll, histotitle=ll+' ;mc_{1}; mn_{1}')
                
                
                
    padMaster = r.TH2F('pm', ';mc_{1}; mn_{1}', 50, 105.0, 300.0, 50, 0.0, 120.0)
    
def convertDsidDict2points(dsidDict={'177501': 0.0, }) :
    values = []
    for d, v in dsidDict.iteritems() :
        mc1, mn1 = dsid2mc1mn1(d)
        values.append((mc1, mn1, v))
    return values

def dsid2mc1mn1(dsid) :
    dsInfo="""
    177501  WH_2Lep_1   130.0  10.0 ;
    177502  WH_2Lep_2   140.0  10.0 ;
    177503  WH_2Lep_3   150.0  0.0  ;
    177504  WH_2Lep_4   152.5  22.5 ;
    177505  WH_2Lep_5   162.5  12.5 ;
    177506  WH_2Lep_6   175.0  0.0  ;
    177507  WH_2Lep_7   165.0  35.0 ;
    177508  WH_2Lep_8   175.0  25.0 ;
    177509  WH_2Lep_9   187.5  12.5 ;
    177510  WH_2Lep_10  200.0  0.0  ;
    177511  WH_2Lep_11  177.5  47.5 ;
    177512  WH_2Lep_12  187.5  37.5 ;
    177513  WH_2Lep_13  200.0  25.0 ;
    177514  WH_2Lep_14  212.5  12.5 ;
    177515  WH_2Lep_15  225.0  0.0  ;
    177516  WH_2Lep_16  190.0  60.0 ;
    177517  WH_2Lep_17  200.0  50.0 ;
    177518  WH_2Lep_18  212.5  37.5 ;
    177519  WH_2Lep_19  225.0  25.0 ;
    177520  WH_2Lep_20  237.5  12.5 ;
    177521  WH_2Lep_21  250.0  0.0  ;
    177522  WH_2Lep_22  202.5  72.5 ;
    177523  WH_2Lep_23  212.5  62.5 ;
    177524  WH_2Lep_24  225.0  50.0 ;
    177525  WH_2Lep_25  237.5  37.5 ;
    177526  WH_2Lep_26  250.0  25.0 ;
    177527  WH_2Lep_27  262.5  12.5 ;
    """
    lines = [l.strip() for l in dsInfo.split(';') if l.strip()]
    lines = [l for l in lines if len(l) and dsid in l]
    if len(lines) :
        tokens = lines[0].strip().split()
        mc1, mn1 = tokens[2], tokens[3]
        return (mc1, mn1)
    else :
        print dsid,' not found'

znValuesAnyes = {
    'mm' : {
        'eq1j' : {
            '177501' : 1.35905
            ,'177502' : 1.23629
            ,'177503' : 0.875974
            ,'177504' : 0.952235
            ,'177505' : 0.689296
            ,'177506' : 0.342558
            ,'177507' : 0.586926
            ,'177508' : 0.527277
            ,'177509' : 0.408544
            ,'177510' : 0.367983
            },
        'ge2j' : {
            '177501' : 1.48882
            ,'177502' : 0.915302
            ,'177503' : 1.00018
            ,'177504' : 0.760198
            ,'177505' : 0.605688
            ,'177506' : 0.391331
            ,'177507' : 0.458379
            ,'177508' : 0.295359
            ,'177509' : 0.271602
            ,'177510' : 0.181077
            },
        },
    'em' : {
        'eq1j' : {
             '177501' : 1.26878
            ,'177502' : 0.800345
            ,'177503' : 0.572663
            ,'177504' : 0.50385
            ,'177505' : 0.436316
            ,'177506' : 0.360832
            ,'177507' : 0.287996
            ,'177508' : 0.478104
            ,'177509' : 0.287256
            ,'177510' : 0.213869
        },
        'ge2j' : {
             '177501' : 1.68595
            ,'177502' : 1.04536
            ,'177503' : 0.954081
            ,'177504' : 0.847949
            ,'177505' : 0.767778
            ,'177506' : 0.599028
            ,'177507' : 0.525896
            ,'177508' : 0.485501
            ,'177509' : 0.460561
            ,'177510' : 0.326818
            }
        }
}

znValuesDavide = {'ee': 
                  {'ge2j':
                  {'177509': 0.07, '177502': 0.34, '177503': -0.2,
                  '177501': 0.12, '177506': -0.08, '177507': 0.0,
                  '177504': 0.12, '177505': -0.02, '177525': -0.19,
                  '177508': 0.04, '177527': -0.23, '177520': -0.15,
                  '177521': -0.13, '177522': -0.01, '177523': -0.14,
                  '177511': -0.0, '177510': -0.01, '177513': 0.04,
                  '177512': -0.08, '177515': -0.14, '177514': -0.14,
                  '177517': -0.18, '177519': -0.09, '177518': -0.22,
                  '177526': -0.32},
                  'eq1j':
                  {'177509': 0.02, '177502': 0.25, '177503': 0.19,
                  '177501': 0.18, '177506': 0.08, '177507': -0.04,
                  '177504': 0.07, '177505': 0.06, '177525': -0.09,
                  '177508': 0.01, '177527': -0.12, '177520': -0.12,
                  '177521': -0.1, '177522': -0.02, '177523': -0.09,
                  '177511': -0.03, '177510': -0.0, '177513': -0.05,
                  '177512': -0.09, '177515': -0.07, '177514': -0.04,
                  '177517': -0.06, '177519': -0.05, '177518': -0.1,
                  '177526': -0.11}},  
                'em':
                {'ge2j':
                {'177509': 0.43, '177502': 0.93, '177503': 0.87,
                '177501': 1.58, '177506': 0.56, '177507': 0.47,
                '177504': 0.74, '177505': 0.75, '177525': 0.06,
                '177508': 0.44, '177527': -0.01, '177520': 0.15,
                '177521': 0.02, '177522': 0.11, '177523': 0.21,
                '177511': 0.51, '177510': 0.32, '177513': 0.11,
                '177512': 0.31, '177515': 0.09, '177514': 0.24,
                '177517': 0.3, '177519': 0.08, '177518': 0.25,
                '177526': 0.13}, 
                'eq1j':
                {'177509': 0.27, '177502': 0.76, '177503': 0.54,
                '177501': 1.21, '177506': 0.34, '177507': 0.27,
                '177504': 0.48, '177505': 0.41, '177525': 0.06,
                '177508': 0.45, '177527': -0.07, '177520': -0.0,
                '177521': -0.03, '177522': 0.12, '177523': 0.06,
                '177511': 0.22, '177510': 0.2, '177513': 0.22,
                '177512': 0.22, '177515': 0.05, '177514': 0.04,
                '177517': 0.16, '177519': 0.05, '177518': 0.13,
                '177526': -0.05}},
                'mm':
                {'ge2j':
                {'177509': 0.3, '177502': 0.97, '177503': 0.96,
                '177501': 1.53, '177506': 0.41, '177507': 0.47,
                '177504': 0.85, '177505': 0.64, '177525': -0.01,
                '177508': 0.34, '177527': -0.1, '177520': 0.07,
                '177521': -0.03, '177522': 0.2, '177523': 0.05,
                '177511': 0.45, '177510': 0.19, '177513': 0.19,
                '177512': 0.3, '177515': 0.06, '177514': 0.09,
                '177517': 0.21, '177519': 0.07, '177518': 0.18,
                '177526': -0.06},
                'eq1j':
                {'177509': 0.43, '177502': 1.3, '177503': 0.92,
                '177501': 1.43, '177506': 0.36, '177507': 0.62,
                '177504': 1.0, '177505': 0.73, '177525': 0.1,
                '177508': 0.56, '177527': -0.02, '177520': 0.13,
                '177521': 0.01, '177522': 0.26, '177523': 0.19,
                '177511': 0.55, '177510': 0.39, '177513': 0.34,
                '177512': 0.4, '177515': 0.16, '177514': 0.26,
                '177517': 0.23, '177519': 0.2, '177518': 0.28,
                '177526': 0.07}}
                }


def plotZns(points=[(0.0, 0.0, 1.0),], histoname='', histotitle='', alsoNegative=False) :
    pm = r.TH2F(histoname if histoname else 'pm', 
                histotitle if histotitle else ';mc_{1}; mn_{1}',
                50, 105.0, 300.0, 50, 0.0, 120.0)
    gr = r.TGraph2D(1)
    gr.SetTitle('g_'+histoname)
    gr.SetMarkerStyle(r.kFullSquare)
    gr.SetMarkerSize(2*gr.GetMarkerSize())
    for i, (x, y, zn) in enumerate(points) :
        firstPoint = gr.GetN()==1 and gr.GetZ()[0]==0.0
        if zn>0.0 or alsoNegative : gr.SetPoint(0 if firstPoint else gr.GetN(),
                                                float(x), float(y), zn)
    c = r.TCanvas('c_zn_'+histoname, histoname, 800, 600)
    c.cd()
    pm.SetStats(0)
    pm.Draw('axis')
    if gr.GetN() : gr.Draw('colz same') #'pcol')
    tex = r.TLatex(0.0, 0.0, '')
    tex.SetTextFont(pm.GetTitleFont())
    tex.SetTextSize(0.75*tex.GetTextSize())
    for x, y, zn in points : tex.DrawLatex(float(x), float(y), "%.2f"%(zn if zn>0.0 or alsoNegative else 0.0))
    c.Update()
    extensions      = ['png']
    for ext in extensions : c.SaveAs(c.GetName()+'.'+ext)

def combineZn(znDictPerChannel={}, alsoNegative=False) :
    """combine the zn from multiple channels by summing them in
    quadrature; return a dict with the same keys (i.e. dsids)
    """
    sqrt = math.sqrt
    allDsets = list(set([k for c,v in znDictPerChannel.iteritems() for k in v.keys()]))
    return dict([(d, sqrt(sum([znThisChan[d]*znThisChan[d]
                               for c, znThisChan in znDictPerChannel.iteritems() if d in znThisChan])))
                  for d in allDsets])


if __name__=='__main__' :
    plot_zn_after_opt()

   
