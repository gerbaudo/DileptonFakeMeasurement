#!/bin/env python

# pick objects from several files, and write them to one output file


import ROOT as r
r.gROOT.SetBatch(1)

input0 = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04/SusyMatrixMethod/data/FinalFakeHist_Apr_10.root'
basedir = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04_fake_dev/SusyTest0/run'
input1 = basedir+'/'+'out/fake_matrix_May_10/fake_matrices_el.root'
input2 = basedir+'/'+'out/fake_matrix_May_10/fake_matrices_mu.root'
input3 = basedir+'/'+'out/fake_matrix_May_10/fake_matrices_el_anygroup.root'
input4 = basedir+'/'+'out/fake_matrix_May_10/fake_matrices_mu_anygroup.root'
input5 = basedir+'/'+'out/fake_matrix_syst_shift_May_10/fake_matrices_el.root'
input6 = basedir+'/'+'out/fake_matrix_syst_shift_May_10/fake_matrices_mu.root'

outfileName = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04/SusyMatrixMethod/data/FinalFakeHist_May_20.root'

objMap = { input0 : [('mu_real_eff_CR_SSInc1j', None) # none = keep the orig name
                     ,('mu_real_eff2d_CR_SSInc1j', None)
                     ,('el_real_eff_CR_SSInc1j', None)
                     ,('el_real_eff2d_CR_SSInc1j', None)
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
                     ,('el_fake_rate_CR_SSInc1j', None) # put it there, but we shouldn't use the 1d param anymore...
                     ,('mu_fake_rate_CR_SSInc1j', None)
                     ,('el_real_eff2d_CR_SRWHnoMlj', None)
                     ,('el_fake_rate2d_CR_SRWHnoMlj', None)
                     ,('mu_real_eff2d_CR_SRWHnoMlj', None)
                     ,('mu_fake_rate2d_CR_SRWHnoMlj', None)
                     ],
           input1 : [('el_fake_pt_eta_ssinc1j', 'el_fake_rate2d_CR_SSInc1j')
                     ],
           input2 : [('mu_fake_pt_eta_ssinc1j', 'mu_fake_rate2d_CR_SSInc1j')
                     ],
           input3 : [('el_fake_pt_eta_ssinc1j', 'el_fake_rate2d_CR_SSInc1j_frac_do')],
           input4 : [('mu_fake_pt_eta_ssinc1j', 'mu_fake_rate2d_CR_SSInc1j_frac_do')],
           input5 : [('el_fake_pt_eta_ssinc1j', 'el_fake_rate2d_CR_SSInc1j_frac_up')],
           input6 : [('mu_fake_pt_eta_ssinc1j', 'mu_fake_rate2d_CR_SSInc1j_frac_up')],
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
        


