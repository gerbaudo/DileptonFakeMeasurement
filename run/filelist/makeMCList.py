# Make lists of input root files

import subprocess
import datasets

validModes = ['data', 'mc12', 'susy']
mode = 'mc12'
assert mode in validModes,"Invalid mode %s (should be one of %s)" % (mode, str(validModes))

# Directory where files are
basedir = {'data' : '/gdata/atlas/ucintprod/SusyNt/data12_n0115/', # data
           'mc12' : '/gdata/atlas/ucintprod/SusyNt/mc12_n0115/',   # mc backgrounds
           'susy' : '/gdata/atlas/ucintprod/SusyNt/susy_n0115/',   # mc signals
           }
# Tag for production
tags = []
tags.append("n0115")
#tags.append("n0111")
wantedDsets = datasets.wantedDsets

###############################################################################################
#                           Don't need to edit below here!!!                                  #
###############################################################################################
dlist = []
dir = basedir[mode]
for tag in tags :
    print tag
    ls = subprocess.Popen(["ls " + dir + " | grep " + tag + " | grep user"],shell=True,stdout=subprocess.PIPE)
    dlist  = dlist + [l for l in [ll.lstrip().rstrip()
                                  for ll in (ls.stdout.read()).split("\n")]
                      if l] # skip empty lines

def contains(dataset, name):
    if (name + ".") in dataset:
        return True
    return False

def makeFile(dataset, name):
    ls = subprocess.Popen(["ls " + dir + dataset + "/* > " + name + ".txt"],shell=True)
    ls.wait()

wanted = wantedDsets[mode]
for ds in dlist:
    print ds
    for name in wanted:
        if(contains(ds,name) and not ("_a" in ds)):
            makeFile(ds,name)


