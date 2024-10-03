import FWCore.ParameterSet.Config as cms

muonsFromScoutingMuons = cms.EDProducer("Run3ScoutingMuonToRecoMuonProducer",

    # Muon collections
    srcMuons = cms.InputTag("hltScoutingMuonPackerVtx")

)
