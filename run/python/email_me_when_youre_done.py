#!/bin/env python

# email me when all my jobs on the queue are done
#
# davide.gerbaudo@gmail.com
# Oct 2013

import optparse
import os
from time import sleep
from utils import getCommandOutput

def main() :
    parser = optparse.OptionParser()
    parser.add_option('--email',        default='davide.gerbaudo@gmail.com')
    parser.add_option('--message',      default='all jobs done, bam!')
    parser.add_option('--polleverysec', default=60*5, type='int')
    parser.add_option('--subject',      default='All jobs done')
    parser.add_option('--user',         default=os.environ['USER'])
    (opts, args) = parser.parse_args() 
    receiver  = opts.email
    username  = opts.user
    message   = opts.message
    subject   = opts.subject
    everyNsec = opts.polleverysec
    nChecksDone = 0
    countJobsCmd="qstat -u %(user)s | grep %(user)s | wc -l"% {'user':username}
    nJobsRunning = int(getCommandOutput(countJobsCmd)['stdout'])    
    while nJobsRunning > 0 :
        sleep(everyNsec)
        nJobsRunning = int(getCommandOutput(countJobsCmd)['stdout'])
        nChecksDone += 1
    message = message+("\n(after %(nchk)d checks every %(nsec)d sec)"
                       %{'nchk' : nChecksDone, 'nsec' : everyNsec})
    mailCmd = ("echo \"%(msg)s\" | mail -s \"%(sbj)s\" %(dest)s"
               % {'dest' : receiver, 'msg' : message, 'sbj':subject})
    getCommandOutput(mailCmd)

if __name__=='__main__' :
    main()
