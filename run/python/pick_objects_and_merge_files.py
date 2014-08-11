#!/bin/env python

# pick objects from several files, and write them to one output file


import ROOT as r
r.gROOT.SetBatch(1)

input0 = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04/SusyMatrixMethod/data/FinalFakeHist_Apr_10.root'
basedir = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04_fake_dev/SusyTest0/run'

input1 = basedir+'/'+'out/fake/weigtedmatrix_Jul_26/fake_matrices_el.root'
input2 = basedir+'/'+'out/fake/weigtedmatrix_Jul_26/fake_matrices_mu.root'

outfileName = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04/SusyMatrixMethod/data/FinalFakeHist_Jul_26.root'
#('orig_name', 'new_name')  # none = keep the orig name
objMap = { input0 : [('mu_real_eff_CR_SSInc1j',    'mu_real_eff_emuInc')
                     ,('mu_real_eff2d_CR_SSInc1j', 'mu_real_eff2d_emuInc')
                     ,('el_real_eff_CR_SSInc1j',   'el_real_eff_emuInc')
                     ,('el_real_eff2d_CR_SSInc1j', 'el_real_eff2d_emuInc')
                     ,('el_real_up', None)
                     ,('el_real_down', None)
                     ,('mu_real_up', None)
                     ,('mu_real_down', None)
                     ,('el_HFLFerr', None)
                     ,('mu_HFLFerr', None)
                     ,('el_datamc', None)
                     ,('mu_datamc', None)
                     ,('el_region', None)
                     ,('mu_region', None)
                     ,('el_eta_sys', None)
                     ,('mu_eta_sys', None)
                     ],
           input1 : [('el_fake_pt_emu',      'el_fake_rate_emuInc')
                     ,('el_fake_pt_eta_emu', 'el_fake_rate2d_emuInc')
                     ],
           input2 : [('mu_fake_pt_emu',      'mu_fake_rate_emuInc')
                     ,('mu_fake_pt_eta_emu', 'mu_fake_rate2d_emuInc')
                     ],
           }
inputFiles = dict((fn, r.TFile.Open(fn, 'read')) for fn in objMap.keys())
outputFile = r.TFile.Open(outfileName, 'recreate')
outputFile.cd()

for fn, objList in objMap.iteritems():
    inFile = inputFiles[fn]
    for oldName, newName in objList:
        newName = newName if newName else oldName
        obj = inFile.Get(oldName)
        if hasattr(obj,'SetName') : obj.SetName(newName)
        obj.Write(newName)
outputFile.Close()
print "objects written to %s"%outfileName
        


