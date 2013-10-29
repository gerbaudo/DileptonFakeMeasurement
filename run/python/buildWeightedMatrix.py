#!/bin/env python

from math import sqrt
import optparse
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options

usage="""
Example usage:
%prog \\
 --tag ${TAG} \\
 --input_dir out/fakerate/merged/data_${TAG}.root \\
 --output_file out/fakerate/merged/FinalFakeHist_${TAG}.root \\
 --output_plot out/fakerate/merged/FinalFakeHist_plots_${TAG} \\
 >& log/fakerate/FinalFakeHist_${TAG}.log
"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-t', '--tag')
    parser.add_option('-i', '--input_dir')
    parser.add_option('-o', '--output_file')
    parser.add_option('-p', '--output_plot')
    parser.add_option('-v','--verbose', action='store_true', default=False)
    (opts, args) = parser.parse_args()
    requiredOptions = ['tag', 'input_dir', 'output_file', 'output_plot']
    otherOptions = ['verbose']
    allOptions = requiredOptions + otherOptions
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions) : parser.error('Missing required option')
    tag = opts.tag
    inputDirname  = opts.input_dir
    outputFname   = opts.output_file
    outputPlotDir = opts.output_plot
    verbose       = opts.verbose
    if verbose : print ('\nUsing the following options:\n'
                        +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))

    inputFiles = getInputFiles(inputDirname, tag)
    assert all(f for f in inputFiles.values()), ("missing inputs: \n%s"
                                                 %'\n'.join(["%s : %s"%kv
                                                             for kv in inputFiles.iteritems()]))
    

def samples() : return ['allBkg', 'ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def getInputFiles(inputDirname, tag) :
    inDir = inputDirname
    tag = tag if tag.startswith('_') else '_'+tag
    return dict(zip(samples(), [r.TFile.Open(inDir+'/'+s+tag+'.root') for s in samples()]))
    
#   FinalNewFake format;
#   format.setTag(tag);
#   format.setInputDir(inputdir);
#   format.setOuputFilename(outputfile);
#   format.setOuputPlotdir(plotdir);
#   format.setDebug(dbg);
#   format.initIoFiles();
#   format.buildRates();

if __name__=='__main__' :
    main()
