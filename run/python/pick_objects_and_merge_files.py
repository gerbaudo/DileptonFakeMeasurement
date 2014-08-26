#!/bin/env python

# pick objects from several files, and write them to one output file
#
# davide.gerbaudo@gmail.com
# 2014

import optparse
import os
import pprint

import ROOT as r
r.gROOT.SetBatch(1)


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input-dict', help='input dictionary, specifying what to pick from where')
    parser.add_option('-o', '--output-file', help='dir for plots')
    parser.add_option('-p', '--print-default-dict', action='store_true', help='print the example input-dict.'
                      "\nFormat: [('orig_name', 'new_name'), ('name1', None)]  # None = keep the orig name")
    (options, args) = parser.parse_args()
    if not options.output_file : parser.error('specify an output file')
    if options.print_default_dict:
        print 'Default input dictionary:'
        pprint.pprint(get_default_input_dict())
        return
    input_dict = eval(open(options.input_dict).read()) if options.input_dict else get_default_input_dict()
    inputFiles = dict((fn, r.TFile.Open(fn, 'read')) for fn in input_dict.keys())
    if not all(v for k,v in inputFiles.iteritems()):
        parser.error('Invalid input files:\n'+pprint.pformat(inputFiles))
    out = r.TFile.Open(options.output_file, 'recreate')
    out.cd()

    for fn, objList in input_dict.iteritems():
        inFile = inputFiles[fn]
        for oldName, newName in objList:
            newName = newName if newName else oldName
            obj = inFile.Get(oldName)
            if not obj:
                print "missing '{0}' from '{1}'".format(oldName, inFile.GetName())
                continue
            if hasattr(obj,'SetName') : obj.SetName(newName)
            obj.Write(newName)
    out.Close()
    print "objects written to %s"%out.GetName()
        

def get_default_input_dict():
    input0 = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04/SusyMatrixMethod/data/FinalFakeHist_Apr_10.root'
    basedir = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04_fake_dev/SusyTest0/run'

    input1 = basedir+'/'+'out/fake/weigtedmatrix_Jul_26/fake_matrices_el.root'
    input2 = basedir+'/'+'out/fake/weigtedmatrix_Jul_26/fake_matrices_mu.root'

    #('orig_name', 'new_name')  # none = keep the orig name
    input_dict = { input0 : [('mu_real_eff_CR_SSInc1j',    'mu_real_eff_emuInc')
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
    return input_dict


if __name__=='__main__':
    main()
