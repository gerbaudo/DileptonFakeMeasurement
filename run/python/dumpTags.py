#!/bin/env python

# Script to loop over the external packages and save their tags to a file.
#
# When a package was checked out without a tag, the revision is saved
# instead. Make sure you have a valid kerberos ticket to avoid having
# to type your password N times. The script extracts the info from two
# lines of the 'svn info' output.
#
# URL: svn+ssh://svn.cern.ch/reps/atlasusr/morgens/TauCorrections/tags/TauCorrections-00-00-07
# Revision: 147016
#
# davide.gerbaudo@gmail.com
# 2013-06-05

import os
import re
import subprocess

path     = '../../'          # assume we are running from SusyTest0/run
outFname = '../doc/tags.txt' # where the tags will be written
def getCommandOutput(command):
    "lifted from supy (https://github.com/elaird/supy/blob/master/utils/io.py)"
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout,stderr = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}
def extractPkg(line) :
    match = re.search('.*/(\S+)/((tags)|(trunk))/?', line)
    return match.group(1) if match else None
def extractRev(line) :
    match = re.search('Revision: (\d+)', line)
    return match.group(1) if match else None
def extractTag(line) :
    match = re.search('.*tags/(\S+)', line)
    return match.group(1) if match else None
def getTagOrRev(dir) :
    "tag if available, otherwise rev"
    cmd = 'svn info '+dir
    res = getCommandOutput(cmd)
    lines = res['stdout'].split('\n')
    urlLine = next((l for l in lines if l.startswith('URL:')),      None)
    revLine = next((l for l in lines if l.startswith('Revision:')), None)
    if not urlLine or not revLine : return
    pkg = extractPkg(urlLine)
    tag = extractTag(urlLine)
    rev = extractRev(revLine)
    assert tag or rev, "cannot extract info for %s"%dir
    return "%s @ %s"%(pkg, tag if tag else rev)


dirs = sorted([path+'/'+p for p in os.listdir(path)])
dirs = filter(os.path.isdir, dirs)
dirs = filter(lambda x : 'SusyTest0' not in x, dirs) # we use git!
labels = []
for d in dirs :
    label = getTagOrRev(d)
    # fallback for the pkgs that appear as 'not a working copy'
    if not label : label = getTagOrRev(d+'/cmt')
    labels.append(label)

with open(outFname, 'w') as outfile :
    outfile.write('\n'.join(labels)+'\n')



