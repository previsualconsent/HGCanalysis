#!/usr/bin/env python

import os,sys
import optparse
import commands


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                                        , default='1nh')
parser.add_option('-t', '--tag'        ,    dest='tag'                , help='tag'                                                , default='MinBias_v4')
parser.add_option('-f', '--ffile'      ,    dest='ffile'              , help='first file'                                         , default=0,  type=int)
parser.add_option('-s', '--step'       ,    dest='step'               , help='step'                                               , default=-1,  type=int)
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'                                   , default='/store/cmst3/group/hgcal/CMSSW/Ntuples/')
parser.add_option('-c', '--cfg'        ,    dest='cfg'                , help='cfg file'                                           , default='test/runHGCSimHitsAnalyzer_cfg.py')
(opt, args) = parser.parse_args()


#prepare output
os.system('cmsMkdir %s'%opt.output)
cmsswBase=os.environ['CMSSW_BASE']
jobsDir=cmsswBase+'/src/FARM'
os.system('mkdir -p %s'%jobsDir)

#create a wrapper for standalone cmssw job
scriptFile = open('%s/runJob_%s_%d.sh'%(jobsDir,opt.tag,opt.ffile), 'w')
scriptFile.write('#!/bin/bash\n')
scriptFile.write('cd %s/src\n'%cmsswBase)
scriptFile.write('eval `scram r -sh`\n')
scriptFile.write('cd %s\n'%jobsDir)
scriptFile.write('cmsRun %s/src/UserCode/HGCanalysis/%s %s %d %d\n'%(cmsswBase,opt.cfg,opt.tag,opt.ffile,opt.step))
scriptFile.write('cmsStage -f /tmp/psilva/%s_SimHits_%d.root %s\n'%(opt.tag,opt.ffile,opt.output))
scriptFile.write('rm /tmp/psilva/%s_SimHits_%d.root\n'%(opt.tag,opt.ffile))
scriptFile.write('echo "All done for this job\n')
scriptFile.close()
os.system('chmod u+rwx %s/runJob_%s_%d.sh'%(jobsDir,opt.tag,opt.ffile))

#submit it to the batch or run it locally
if opt.queue=='':
    os.system('%s/runJob_%s_%d.sh'%(jobsDir,opt.tag,opt.ffile))
else:
    os.system("bsub -q %s \'%s/runJob_%s_%d.sh\'"%(opt.queue,jobsDir,opt.tag,opt.ffile))
    
