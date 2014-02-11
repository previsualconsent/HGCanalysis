# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleElectronPt35_cfi --conditions auto:startup -s GEN,SIM --datatier GEN-SIM -n 10 --relval 9000,100 --eventcontent RAWSIM --geometry Extended2023HGCal --conditions auto:upgradePLS3 --fileout file:singleele_pt35.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

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
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt35_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/tmp/Events_XXX_SEED_XXX.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
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

process.generator = cms.EDProducer("Pythia6PtYDistGun",
				   PGunParameters = cms.PSet( MaxY = cms.double(5),
							      MinY = cms.double(-5),
							      YBinning = cms.int32(500),
							      MinPt = cms.double(5.0),
							      MaxPt = cms.double(500.0),
							      PtBinning = cms.int32(100000),
							      MinPhi = cms.double(-3.1415),
							      MaxPhi = cms.double(3.1415),
							      ParticleID = cms.vint32(11),
							      kinematicsFile = cms.FileInPath('UserCode/HGCanalysis/test/particle_gun_pdf.root')
							      ),
				   pythiaPylistVerbosity = cms.untracked.int32(1),
				   firstRun = cms.untracked.uint32(1),
				   Verbosity = cms.untracked.int32(0),
				   psethack = cms.string('electron gun'),
				   pythiaHepMCVerbosity = cms.untracked.bool(False),
				   maxEventsToPrint = cms.untracked.int32(5),
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
								parameterSets = cms.vstring('pythiaUESettings')
								)
				   )




# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

