

from rootUtils import importRoot, buildRatioHistogram
r = importRoot()

def samples() : return ['allBkg', 'ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def fakeProcesses() : return ['ttbar', 'wjets', 'zjets', 'diboson', 'heavyflavor']
def getInputFiles(inputDirname, tag, verbose=False) :
    inDir = inputDirname
    tag = tag if tag.startswith('_') else '_'+tag
    files = dict(zip(samples(), [r.TFile.Open(inDir+'/'+s+tag+'.root') for s in samples()]))
    if verbose : print "getInputFiles('%s'):\n\t%s"%(inputDirname, '\n\t'.join("%s : %s"%(k, f.GetName()) for k, f in files.iteritems()))
    return files
def buildRatio(inputFile=None, histoBaseName='') :
    num, den = inputFile.Get(histoBaseName+'_num'), inputFile.Get(histoBaseName+'_den')
    return buildRatioHistogram(num, den, histoBaseName +'_rat')
