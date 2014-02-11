import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.Geometry.GeometryExtended2023HGCal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.TFileService = cms.Service("TFileService", fileName = cms.string('HGCSimHitsAnalysis.root') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from UserCode.HGCanalysis.storeTools_cff import fillFromStore

process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring('file:particlegun.root')
                            #fileNames=cms.untracked.vstring('file:minbias.root')
                            )
process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/MinBias/v0')
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.hgcSimHitsAnalyzer = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                                            ddViewName     = cms.untracked.string(""),
                                            eeHits         = cms.untracked.string("HGCHitsEE"),
                                            heHits         = cms.untracked.string("HGCHitsHE"),
                                            genSource      = cms.untracked.string("genParticles"),
                                            )

process.p = cms.Path(process.hgcSimHitsAnalyzer)

