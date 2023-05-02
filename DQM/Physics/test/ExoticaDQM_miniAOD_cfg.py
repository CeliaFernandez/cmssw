import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


## VarParsing object

options = VarParsing('analysis')

options.register('globaltag', '113X_mcRun3_2021_realistic_v4', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Set Global Tag')
options.register('name', 'TEST', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Set process name')

options.parseArguments()


## Process

process = cms.Process("ExoticaDQM")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = cms.string(options.globaltag)

process.load("DQM.Physics.ExoticaDQM_miniAOD_cfi")
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")


process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.DQMStore.verbose = cms.untracked.int32(2)
process.dqmSaver.workflow = cms.untracked.string('/Physics/ExoticaMiniAOD/' + options.name)
process.dqmSaver.convention = cms.untracked.string('RelVal')


## output
#Needed to the the plots locally

process.output = cms.OutputModule("PoolOutputModule",
  fileName       = cms.untracked.string('EXOTICA_miniAOD_DQM_' + options.name + '.root'),
  outputCommands = cms.untracked.vstring(
    'drop *_*_*_*',
    'keep *_*_*_ExoticaDQM'
    ),
  splitLevel     = cms.untracked.int32(0),
  dataset = cms.untracked.PSet(
    dataTier   = cms.untracked.string(''),
    filterName = cms.untracked.string('')
  )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
       options.inputFiles       
    )
)

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')


## path definitions
process.p      = cms.Path(
    process.ExoticaDQM_miniAOD

)

process.endjob = cms.Path(
    process.endOfProcess
)

process.fanout = cms.EndPath(
    process.output
)


## schedule definition
process.schedule = cms.Schedule(
    process.p,
    process.endjob,
    process.fanout
)
                          


