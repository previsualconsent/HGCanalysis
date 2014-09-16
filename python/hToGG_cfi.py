import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         crossSection = cms.untracked.double(1),
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(14000.0),
                         PythiaParameters = cms.PSet( processParameters = cms.vstring(
                             'Main:timesAllowErrors = 10000',
                             'ParticleDecays:limitTau0 = on',
                             'ParticleDecays:tauMax = 10',
                             'Tune:ee 3',
                             'Tune:pp 5',
                             'HiggsSM:gg2H = on',
                             '25:onMode = off',      # turn OFF all H decays
                             '25:onIfMatch = 22 22', # turn ON H->gamma gamma
                             ),
                             parameterSets = cms.vstring('processParameters')
                             )
                         )
