#!/usr/bin/env python

import os,sys
import optparse
import commands
import time


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                                        , default='1nh')
parser.add_option('-n', '--njobs'      ,    dest='njobs'              , help='number of jobs'                                     , default=1,  type=int)
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'                                   , default='/store/cmst3/group/hgcal/CMSSW/MinBias')
parser.add_option('-c', '--cfg'        ,    dest='cfg'                , help='cfg file'                                           , default='test/runMinBias_GEN_SIM_cfg.py')
parser.add_option('-r', '--rep'        ,    dest='customReplacements' , help='sed replacements for cfg  key1:val1,key2:val2,...'  , default=None)
(opt, args) = parser.parse_args()


#prepare output
os.system('cmsMkdir %s'%opt.output)
cmsswBase=os.environ['CMSSW_BASE']
jobsDir=cmsswBase+'/src/FARM%s'%(time.time())
os.system('mkdir -p %s'%jobsDir)

def replfunc(match):
    return repldict[match.group(0)]


#loop over the required number of jobs
for n in xrange(0,opt.njobs):

    jobSeed=n+1

    #sed the cfg template 
    inCfg = open(opt.cfg).read()
    outCfg = open('%s/cmssw_%d_cfg.py'%(jobsDir,jobSeed), 'w')
    replacements = {'XXX_SEED_XXX':str(jobSeed)}
    if opt.customReplacements is not None:
        for rep in opt.customReplacements.split(','):
            repKeys=rep.split(':')
            replacements[repKeys[0]]=repKeys[1]
    
    for i in replacements.keys():
        inCfg = inCfg.replace(i, replacements[i])
    outCfg.write(inCfg)
    outCfg.close()
    
    #create a wrapper for standalone cmssw job
    scriptFile = open('%s/runJob_%d.sh'%(jobsDir,jobSeed), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('cd %s/src\n'%cmsswBase)
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd %s\n'%jobsDir)
    scriptFile.write('cmsRun cmssw_%d_cfg.py\n'%jobSeed)
    scriptFile.write('cmsStage -f /tmp/Events_%d.root %s\n'%(jobSeed,opt.output))
    scriptFile.write('rm /tmp/Events_%d.root\n'%jobSeed)
    scriptFile.write('echo "All done for job %d" \n'%jobSeed)
    scriptFile.close()
    os.system('chmod u+rwx %s/runJob_%d.sh'%(jobsDir,jobSeed))

    #submit it to the batch or run it locally
    if opt.queue=='':
        os.system('%s/runJob_%d.sh'%(jobsDir.jobSeed))
    else:
        os.system("bsub -q %s -J HGCSIM%d \'%s/runJob_%d.sh\'"%(opt.queue,jobSeed,jobsDir,jobSeed))
    
