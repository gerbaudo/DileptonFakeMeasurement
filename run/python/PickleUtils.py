#!/bin/env python

# Utility functions to read/write stuff to pickle files
#
# davide.gerbaudo@gmail.com
# Jan 2013

import pickle

def dumpToPickle(filename='', obj=None) :
    output = open(filename, 'wb')
    pickle.dump(obj, output)
    output.close()

