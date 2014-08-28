import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCGeometryAnalyzer",
                          geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive')
                          )
