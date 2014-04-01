#!/bin/env python

# extract from the logs the printouts of the suspicious events in the fake control region
#
# Usage:
# grep_suspicious.py log/susyplot/*period*${TAG}.log
#
# davide.gerbaudo@gmail.com
# 2014-04-01

import sys

def main():
    for filename in sys.argv[1:]:
        with open(filename) as f:
            extractText(f)

def extractText(file):
    """Extract a block that looks like this:
    fake suspicious event run 210302 event 73750293 em 2j
    Run 210302 Event 73750293 Stream Muons w 1
    Mu : q  1 pt  46.31 eta  1.98 phi  2.69 cb 1 type 0 origin 0
    El : q  1 pt  10.47 eta  0.14 phi  1.66 tight 1 type 0 origin 0
    Met : pt  35.48 phi 0.65
    Jet : pt  59.82 eta  0.72 phi -1.12 jvf 0.97 mv1  0.05
    Jet : pt  27.30 eta  0.67 phi  2.44 jvf 0.47 mv1  0.06
    Jet : pt  21.93 eta -1.83 phi  1.34 jvf 0.07 mv1  0.27
    """
    inBlock = False
    def startBlock(l) : return l.startswith('fake suspicious event')
    def stillInBlock(l) : return any(l.startswith(t) for t in ['Run ', 'Mu :', 'El :', 'Met :', 'Jet :'])
    buffer = []
    for l in file.readlines():
        if not inBlock and startBlock(l) :
            print ''.join(buffer)
            inBlock = True
            buffer = [l]
        elif inBlock and stillInBlock(l) :
            buffer.append(l)
        else:
            inBlock = False
    if buffer:
        print ''.join(buffer)
        
    

if __name__=='__main__':
    main()
