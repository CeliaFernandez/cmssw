import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
ExoticaDQM_miniAOD = DQMEDAnalyzer(
    "ExoticaDQM_miniAOD",

    #Trigger Results
    TriggerResults           = cms.InputTag('TriggerResults','','HLT'),
    HltPaths                 = cms.vstring("HLT_Mu","HLT_Ele","HLT_Photon","HLT_PFHT","HLT_HT","HLT_PFMET","HLT_MET","HLT_"),

    #Physics objects
    vertexCollection         = cms.InputTag('offlineSlimmedPrimaryVertices'),
    electronCollection       = cms.InputTag("slimmedElectrons"),

    muonCollection           = cms.InputTag("slimmedMuons"),

    photonCollection         = cms.InputTag("slimmedPhotons"),

    DiJetPFJetCollection     = cms.VInputTag('slimmedJets','slimmedJetsAK8'),

    MonoJetPFJetCollection   = cms.VInputTag('slimmedJets','slimmedJetsAK8'),

    pfMETCollection          = cms.InputTag("slimmedMETs"),

    trackCollection          = cms.InputTag("isolatedTracks"),

    displacedSAMuonCollection  = cms.InputTag("displacedStandAloneMuons"),

    # MC truth
    genParticleCollection    = cms.InputTag("prunedGenParticles"),
    
    #Cuts
    # DiJet
    dijet_PFJet1_pt_cut       = cms.double(30.0),
    dijet_PFJet2_pt_cut       = cms.double(30.0),
    # DiMuon
    dimuon_Muon1_pt_cut      = cms.double(30.0),
    dimuon_Muon2_pt_cut      = cms.double(30.0),
    # DiElectron
    dielectron_Electron1_pt_cut = cms.double(30.0),
    dielectron_Electron2_pt_cut = cms.double(30.0),
    # DiPhoton
    diphoton_Photon1_pt_cut   = cms.double(20.0),
    diphoton_Photon2_pt_cut   = cms.double(20.0),
    # MonoMuon
    monomuon_Muon_pt_cut      = cms.double(30.0),
    monomuon_Muon_met_cut     = cms.double(60.0),
    # MonoElectron
    monoelectron_Electron_pt_cut  = cms.double(30.0),
    monoelectron_Electron_met_cut = cms.double(60.0),
    # Monojet
    monojet_PFJet_pt_cut      = cms.double(30.0),
    monojet_PFJet_met_cut     = cms.double(60.0),
    # Displaced lepton or jet
    dispFermion_eta_cut = cms.double(2.4),
    dispFermion_pt_cut  = cms.double(1.0)
    

)
