import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

ffile=0
step=-1
preFix='MinBias_v4'

#configure from command line
import sys
if(len(sys.argv)>2):
    preFix=sys.argv[2]
    if(len(sys.argv)>3):
        if(sys.argv[3].isdigit()) : ffile=int(sys.argv[3])
    if(len(sys.argv)>4):
        if(sys.argv[4].isdigit()) : step=int(sys.argv[4])

process.TFileService = cms.Service("TFileService", fileName = cms.string('/tmp/psilva/%s_SimHits_%d.root'%(preFix,ffile)))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

from UserCode.HGCanalysis.storeTools_cff import fillFromStore

process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring()
                            )
process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)
#process.source.fileNames=cms.untracked.vstring('file:///tmp/psilva/Events_485.root')

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

print process.source.fileNames

process.hgcSimHitsAnalyzer = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                                            ddViewName     = cms.untracked.string(""),
                                            hitCollections = cms.untracked.vstring('HGCHitsEE',  'HGCHitsHEfront'),#      'HGCHitsHEback'           ),
                                            sdTags         = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive'),#  'HGCalHEScintillatorSensitive'),
                                            genSource      = cms.untracked.string("genParticles"),
                                            )

process.p = cms.Path(process.hgcSimHitsAnalyzer)

