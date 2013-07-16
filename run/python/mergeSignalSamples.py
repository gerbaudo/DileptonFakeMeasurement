#!/bin/env python

# merge signal samples to increase the signal stats
#
# davide.gerbaudo@gmail.com
# Jun 2013

import collections
import glob
import os
import optparse
import re
import sys
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(1)

from NavUtils import getAllHistoNames
from SampleUtils import guessSampleFromFilename, ModeAWhDbPar, ModeAWhDbReqid, ModeAWhDbMergedFake2Lreqid

class SamplesMerger :
    def __init__(self) :
        self.filenamesByReqid = collections.defaultdict(list)
        self.masspointByReqid = collections.defaultdict(list)
        self.reqidProvider = ModeAWhDbMergedFake2Lreqid()
        self.parDb = ModeAWhDbPar()
        self.reqDb = ModeAWhDbReqid()
        self.overwrite = False
        self.verbose = False
        self.regexHist = '.*'
    def idTagFromFilename(self, fname) :
        'parse something like wA_noslep_WH_2Lep_9_May16_n0139.AnaHists'
        match = re.search('wA_noslep_WH_2Lep_(?P<id>\d+?)_(?P<tag>.*?).AnaHist', fname)# nongreedy
        return match.group('id'), match.group('tag')
    def addFile(self, filename) :
        sample = guessSampleFromFilename(filename)
        mc1, mn1 = self.parDb.mc1Mn1ByReqid(self.reqDb.reqidBySample(sample))
        fakeReqid = self.reqidProvider.reqidByMc1Mn1(mc1, mn1)
        self.filenamesByReqid[fakeReqid].append(filename)
        self.masspointByReqid[fakeReqid].append((mc1, mn1))
    def mergeAndWrite(self, outdir) :
        """write merged histos, unless !overwrite, in which case just
        compute merged coords and fnames"""
        def mergeFiles(targetFile, filenames, regex) :
            infiles = [r.TFile.Open(f) for f in filenames]
            histonames = [h for h in getAllHistoNames(infiles[0]) if re.search(regex, h)]
            scale = 1.0/len(infiles)
            outFile = r.TFile.Open(targetFile, 'recreate')
            outFile.cd()
            for hn in histonames :
                h = infiles[0].Get(hn).Clone()
                h.Scale(scale)
                for inf in infiles[1:] : h.Add(inf.Get(hn), scale)
                h.Write()
            outFile.Close()
        self.mergedPointsByReqid = {}
        self.outFnameByReqid     = {}
        reqids = self.filenamesByReqid.keys()
        for rid in reqids :
            fnames = self.filenamesByReqid[rid]
            points = self.masspointByReqid[rid]
            npoints = len(points)
            assert npoints,"cannot merge 0 samples"
            def avgPoint(points) :
                np = float(len(points))
                return (sum([x for x,y in points])/np,
                        sum([y for x,y in points])/np)
            point = avgPoint(points)
            ids, tags = zip(*[self.idTagFromFilename(f) for f in fnames]) # vertical slice
            assert len(set(tags))==1,"cannot merge samples with different tags %s"%str(fnames)
            outfname = outdir+'/'+os.path.basename(fnames[0]).replace(str(ids[0]), str(rid))
            self.outFnameByReqid[rid], self.mergedPointsByReqid[rid] = outfname, point
            if self.verbose : print str(point)+' : '+outfname
            skipMerge = os.path.exists(outfname) and not self.overwrite
            if  skipMerge : continue
            else          : mergeFiles(outfname, fnames, self.regexHist)
    def printMergedFiles(self) :
        for rid in self.outFnameByReqid.keys() :
            print self.mergedPointsByReqid[rid], ' : ', self.outFnameByReqid[rid]

        
                

if __name__=='__main__' :
    usage="""Merge the signal histograms produced by SusyPlot.
    For now, only working for WH_2Lep.
    Examples:
    > ./python/mergeSignalSamples.py inputdir/ outputdir/ tag
    > ./python/mergeSignalSamples.py anaplots/ mergedSignals/ Jun06_n0139 -v -r '.*onebin.*'
    """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-o', '--overwrite', action='store_true', dest='overwrite', default=False,
                      help='overwrite existing files')
    parser.add_option('-r', '--regex-histo', dest='regexhist', default='.*',
                      help='only merge the histos whose name matches regexp (runs faster)')
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False,
                      help='print more details about what is going on')
    (options, args) = parser.parse_args()
    if len(args)!=3 : parser.error('incorrect number of arguments')
    inputDir, outputDir, tag = args[0], args[1], args[2]
    overwrite, regexHist, verbose = options.overwrite, options.regexhist, options.verbose
    
    inputFileNames = glob.glob(inputDir+'/'+'*WH_2Lep*'+tag+'*.root')
    if not os.path.isdir(outputDir) : os.mkdir(outputDir)
    if verbose :
        print 'Options:\n' \
          + '\n'.join(["%s : %s" % (o, eval(o))
                       for o in ['inputDir', 'outputDir', 'tag',]])
        print 'Input files:\n'+'\n'.join(inputFileNames)
    
    merger = SamplesMerger()
    merger.regexHist, merger.verbose = regexHist, verbose
    for f in inputFileNames : merger.addFile(f)
    if overwrite : merger.overwrite = True
    merger.mergeAndWrite(outputDir) 
    if verbose : merger.printMergedFiles()
