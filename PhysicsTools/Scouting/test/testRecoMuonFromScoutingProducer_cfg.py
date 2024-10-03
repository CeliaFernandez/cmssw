import FWCore.ParameterSet.Config as cms
process = cms.Process("SCOUT")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
 

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v4')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    "file:/eos/cms/store/data/Run2024G/ScoutingPFMonitor/RAW/v1/000/385/152/00000/085e5565-2390-4f55-9682-7b597d2b56ec.root"
    )
)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('output.root'))
#                               outputCommands = cms.untracked.vstring( 'drop *',
#                                                                       'keep *_hltFEDSelectorL1_*_*',
#                                                                       'keep *_hltScoutingEgammaPacker_*_*',
#                                                                       'keep *_hltScoutingMuonPackerNoVtx_*_*',
#                                                                       'keep *_hltScoutingMuonPackerVtx_*_*',
#                                                                       'keep *_hltScoutingMuonPacker_*_*',
#                                                                       'keep *_hltScoutingPFPacker_*_*',
#                                                                       'keep *_hltScoutingPrimaryVertexPacker_*_*',
#                                                                       'keep *_hltScoutingTrackPacker_*_*',
#                                                                       'keep *_genParticles_*_*',
#                                                                       'keep edmTriggerResults_*_*_*' ))


process.load("PhysicsTools.Scouting.muonsFromScoutingMuons_cfi")
   
process.p = cms.Path(process.muonsFromScoutingMuons)

process.e = cms.EndPath(process.out)
