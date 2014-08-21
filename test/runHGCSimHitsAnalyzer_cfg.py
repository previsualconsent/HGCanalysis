import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

# configure from command line
# cmsRun test/runHGCSimHitsAnalyzer_cfg.py tag
# where tag can be any sub-directory under /store/cmst3/group/hgcal/CMSSW
#           or any upgrade relval sample (may need tweaking for new releases...)
ffile=0
step=-1
preFix='Single211_SLHC16'
import os,sys
if(len(sys.argv)>2):
    preFix=sys.argv[2]
    if(len(sys.argv)>3):
        if(sys.argv[3].isdigit()) : ffile=int(sys.argv[3])
    if(len(sys.argv)>4):
        if(sys.argv[4].isdigit()) : step=int(sys.argv[4])

#configure the source (list all files in directory within range [ffile,ffile+step[
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring()
                            )
if preFix.find('RelVal')>=0 :
    cmsswVersion=os.environ['CMSSW_VERSION']
    process.source.fileNames=fillFromStore('/store/relval/%s/%s/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/'%(cmsswVersion,preFix),ffile,step)
else :
    process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#run the analyzer
import getpass
whoami=getpass.getuser()
process.TFileService = cms.Service("TFileService", fileName = cms.string('/tmp/%s/%s_SimHits_%d.root'%(whoami,preFix,ffile)))
process.load('UserCode.HGCanalysis.hgcSimHitsAnalyzer_cfi')
process.p = cms.Path(process.analysis)

