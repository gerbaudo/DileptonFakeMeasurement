# Make lists of input root files

import subprocess

# Directory where files are
dir = "/gdata/atlas/mrelich/SusySkims/"

# What files we want, with appropriate names
periods = ("A", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12")
streams = ("Egamma", "Muons")


###############################################################################################
#                           Don't need to edit below here!!!                                  #
###############################################################################################


ls = subprocess.Popen(["ls " + dir + " | grep user"],shell=True,stdout=subprocess.PIPE)
dlist = (ls.stdout.read()).split("\n")
del dlist[-1]

def contains(dataset, name, stream):
    if name in dataset:
        if stream in dataset:
            return True
    return False

def makeFile(dataset, name):
    ls = subprocess.Popen(["ls " + dir + dataset + "/* > " + name + ".txt"],shell=True)
    ls.wait()

for ds in dlist:
    for period in periods:
        per = "period" + period
        for stream in streams:
            if(contains(ds,per,stream)):
                outds = per + "_phs_" + stream           
                makeFile(ds,outds)

    
