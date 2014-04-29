import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#preFix='MinBias_v4'
preFix='SingleMuon_v5'
process.TFileService = cms.Service("TFileService", fileName = cms.string('/data/psilva/%s_SimHits.root'%preFix) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from UserCode.HGCanalysis.storeTools_cff import fillFromStore

process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring()
                            )
process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

print process.source.fileNames

process.hgcSimHitsAnalyzer = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                                            ddViewName     = cms.untracked.string(""),
                                            hitCollections = cms.untracked.vstring('HGCHitsEE',  'HGCHitsHEfront',      'HGCHitsHEback'           ),
                                            sdTags         = cms.untracked.vstring('EESensitive','HESiliconSensitive',  'HEScintillatorSensitive'),
                                            cellSizePars   = cms.untracked.vint32 ( 0,            1,                   1                       ),
                                            genSource      = cms.untracked.string("genParticles"),
                                            )

process.p = cms.Path(process.hgcSimHitsAnalyzer)

