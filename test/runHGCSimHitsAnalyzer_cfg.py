import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.TFileService = cms.Service("TFileService", fileName = cms.string('HGCSimHitsAnalysis.root') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from UserCode.HGCanalysis.storeTools_cff import fillFromStore

process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring('file:particlegun.root')
                            #fileNames=cms.untracked.vstring('file:minbias.root')
                            )
#process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/MinBias_v1')
#process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/SingleMuon_v1')
#process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/SingleElectron_v1')
process.source.fileNames=cms.untracked.vstring('file:/afs/cern.ch/user/p/psilva/work/CMSSW_6_2_0_SLHC8/src/UserCode/HGCanalysis/Events_1.root',
                                               'file:/afs/cern.ch/user/p/psilva/work/CMSSW_6_2_0_SLHC8/src/UserCode/HGCanalysis/Events_2.root',
                                               'file:/afs/cern.ch/user/p/psilva/work/CMSSW_6_2_0_SLHC8/src/UserCode/HGCanalysis/Events_3.root',
                                               'file:/afs/cern.ch/user/p/psilva/work/CMSSW_6_2_0_SLHC8/src/UserCode/HGCanalysis/Events_4.root',
                                               'file:/afs/cern.ch/user/p/psilva/work/CMSSW_6_2_0_SLHC8/src/UserCode/HGCanalysis/Events_1847.root'
                                               )
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

