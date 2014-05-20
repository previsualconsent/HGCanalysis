import FWCore.ParameterSet.Config as cms

process = cms.Process('SIMDIGI')

usePythia8=True

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(250)
)

process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)
process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string(''),
    annotation = cms.untracked.string('runParticleGun_GEN_SIM_cfg nevts:250'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('/tmp/Events_XXX_SEED_XXX.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

if usePythia8:
	process.generator = cms.EDFilter("Pythia8GeneratorFilter",
					 crossSection = cms.untracked.double(1),
					 maxEventsToPrint = cms.untracked.int32(0),
					 pythiaPylistVerbosity = cms.untracked.int32(1),
					 filterEfficiency = cms.untracked.double(1.0),
					 pythiaHepMCVerbosity = cms.untracked.bool(False),
					 comEnergy = cms.double(14000.0),
					 PythiaParameters = cms.PSet( processParameters = cms.vstring( 'Main:timesAllowErrors = 10000',
												       'ParticleDecays:limitTau0 = on',
												       'ParticleDecays:tauMax = 10',
												       'SoftQCD:all = on',
												       'Tune:ee 3',
												       'Tune:pp 5'),
								      parameterSets = cms.vstring('processParameters')
								      )
					 )
else :
	process.generator = cms.EDFilter("Pythia6GeneratorFilter",
					 maxEventsToPrint = cms.untracked.int32(0),
					 pythiaPylistVerbosity = cms.untracked.int32(0),
					 filterEfficiency = cms.untracked.double(1.0),
					 pythiaHepMCVerbosity = cms.untracked.bool(False),
					 comEnergy = cms.double(14000.0),
					 crossSection = cms.untracked.double(1),
					 UseExternalGenerators = cms.untracked.bool(False),
					 PythiaParameters = cms.PSet( pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
												     'MSTJ(22)=2     ! Decay those unstable particles', 
												     'PARJ(71)=10 .  ! for which ctau  10 mm', 
												     'MSTP(33)=0     ! no K factors in hard cross sections', 
												     'MSTP(2)=1      ! which order running alphaS', 
												     'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
												     'MSTP(52)=2     ! work with LHAPDF', 
												     'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
												     'PARP(89)=1800. ! sqrts for which PARP82 is set', 
												     'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
												     'MSTP(95)=6     ! CR (color reconnection parameters)', 
												     'PARP(77)=1.016 ! CR', 
												     'PARP(78)=0.538 ! CR', 
												     'PARP(80)=0.1   ! Prob. colored parton from BBR', 
												     'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
												     'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
												     'PARP(62)=1.025 ! ISR cutoff', 
												     'MSTP(91)=1     ! Gaussian primordial kT', 
												     'PARP(93)=10.0  ! primordial kT-max', 
												     'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
												     'MSTP(82)=4     ! Defines the multi-parton model'),
								      processParameters = cms.vstring('MSEL=0 ! User defined processes',
												      'MSUB(11)=1 ! Min bias process',
												      'MSUB(12)=1 ! Min bias process',
												      'MSUB(13)=1 ! Min bias process',
												      'MSUB(28)=1 ! Min bias process',
												      'MSUB(53)=1 ! Min bias process',
												      'MSUB(68)=1 ! Min bias process',
												      'MSUB(92)=1 ! Min bias process, single diffractive',
												      'MSUB(93)=1 ! Min bias process, single diffractive',
												      'MSUB(94)=1 ! Min bias process, double diffractive',
												      'MSUB(95)=1 ! Min bias process'),
								      parameterSets = cms.vstring('pythiaUESettings', 'processParameters')
								      )
					 )
	
	
# Other statements
from SimGeneral.MixingModule.hgcalDigitizer_cfi import *
process.theDigitizersValid.hgceeDigitizer=hgceeDigitizer
process.theDigitizersValid.hgchebackDigitizer=hgchebackDigitizer
process.theDigitizersValid.hgchefrontDigitizer=hgchefrontDigitizer
process.mix.digitizers=cms.PSet(process.theDigitizersValid)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,
				process.digitisation_step,process.L1simulation_step,process.digi2raw_step,
				process.endjob_step,process.output_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023NoEE

#call to customisation function cust_2023NoEE imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023NoEE(process)

