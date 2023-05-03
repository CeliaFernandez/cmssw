#ifndef ExoticaDQM_miniAOD_H
#define ExoticaDQM_miniAOD_H

#include <memory>

// DQM
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DataKeyTags.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
// Other
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// PAT candidates
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

// ROOT
#include "TLorentzVector.h"

// STDLIB
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>

class ExoticaDQM_miniAOD : public DQMEDAnalyzer {
public:
  ExoticaDQM_miniAOD(const edm::ParameterSet& ps);
  ~ExoticaDQM_miniAOD() override;

protected:
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;

  //Resonances
  virtual void analyzeDiJets(edm::Event const& e);
  virtual void analyzeDiMuons(edm::Event const& e);
  virtual void analyzeDiElectrons(edm::Event const& e);
  virtual void analyzeDiPhotons(edm::Event const& e);

  //Mono Searches
  virtual void analyzeMonoJets(edm::Event const& e);
  virtual void analyzeMonoMuons(edm::Event const& e);
  virtual void analyzeMonoElectrons(edm::Event const& e);

  // Displaced Fermion Searches
  virtual void analyzeDisplacedLeptons(edm::Event const& e, const edm::EventSetup& s);
  virtual void analyzeDisplacedJets(edm::Event const& e, const edm::EventSetup& s);

  // Estimate the momentum vector that a GenParticle would have at its trajectory's point of closest approach to the beam-line.
  virtual GlobalVector getGenParticleTrajectoryAtBeamline(const edm::EventSetup& iSetup, const reco::GenParticle* gen);

private:
  void bookHistograms(DQMStore::IBooker& bei, edm::Run const&, edm::EventSetup const&) override;

  int nLumiSecs_;
  int nEvents_, irun, ievt;

  //Vertex
  edm::EDGetTokenT<reco::VertexCollection> VertexToken_;
  edm::Handle<reco::VertexCollection> VertexCollection_;

  // Electrons
  edm::EDGetTokenT<std::vector<pat::Electron> > ElectronToken_;
  edm::Handle<std::vector<pat::Electron> > ElectronCollection_;

  // Muons
  edm::EDGetTokenT<std::vector<pat::Muon> > MuonToken_;
  edm::Handle<std::vector<pat::Muon> > MuonCollection_;

  // Photons
  edm::EDGetTokenT<std::vector<pat::Photon> > PhotonToken_;
  edm::Handle<std::vector<pat::Photon> > PhotonCollection_;

  // DiJet Jets
  std::vector<edm::EDGetTokenT<std::vector<pat::Jet> > > DiJetPFJetToken_;
  std::vector<edm::InputTag> DiJetPFJetCollection_;
  edm::Handle<std::vector<pat::Jet> > DiJetpfJetCollection_;
  std::vector<pat::Jet> DiJetpfjets;

  // MonoJet Jets
  std::vector<edm::EDGetTokenT<std::vector<pat::Jet> > > MonoJetPFJetToken_;
  std::vector<edm::InputTag> MonoJetPFJetCollection_;
  edm::Handle<std::vector<pat::Jet> > MonoJetpfJetCollection_;
  std::vector<pat::Jet> MonoJetpfjets;

  // MET
  edm::EDGetTokenT<std::vector<pat::MET> > PFMETToken_;
  edm::Handle<std::vector<pat::MET> > pfMETCollection_;

  // Tracks
  edm::EDGetTokenT<std::vector<pat::IsolatedTrack> > TrackToken_;
  edm::Handle<std::vector<pat::IsolatedTrack> > TrackCollection_;

  // Special collections for highly displaced particles
  edm::EDGetTokenT<reco::TrackCollection> MuonDispSAToken_;
  edm::Handle<reco::TrackCollection> MuonDispSACollection_;
  edm::EDGetTokenT<reco::TrackCollection> MuonDispGLToken_;
  edm::Handle<reco::TrackCollection> MuonDispGLCollection_;

  // MC truth
  edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
  edm::Handle<reco::GenParticleCollection> GenCollection_;

  //ES tokens
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;

  ///////////////////////////
  // Parameters
  ///////////////////////////
  // Cuts - MultiJets
  // inputs

  //Varibles Used
  // PFJets
  double PFJetPx[2];
  double PFJetPy[2];
  double PFJetPt[2];
  double PFJetEta[2];
  double PFJetPhi[2];
  double PFJetRapidity[2];
  double PFJetMass[2];
  double PFJetNHEF[2];
  double PFJetCHEF[2];
  double PFJetNEMF[2];
  double PFJetCEMF[2];

  // Muons
  //
  double MuonPx[2];
  double MuonPy[2];
  double MuonPt[2];
  double MuonEta[2];
  double MuonPhi[2];
  double MuonCharge[2];

  // Electrons
  //
  double ElectronPx[2];
  double ElectronPy[2];
  double ElectronPt[2];
  double ElectronEta[2];
  double ElectronPhi[2];
  double ElectronCharge[2];

  // Photon
  //
  double PhotonEnergy[2];
  double PhotonPt[2];
  double PhotonEt[2];
  double PhotonEta[2];
  double PhotonEtaSc[2];
  double PhotonPhi[2];
  double PhotonHoverE[2];
  double PhotonSigmaIetaIeta[2];
  double PhotonTrkSumPtSolidConeDR03[2];
  double PhotonE1x5E5x5[2];
  double PhotonE2x5E5x5[2];

  // MET
  //
  double PFMETEt;
  double PFMETPhi;


  ///////////////////////////
  // Histograms
  ///////////////////////////
  // Histograms - Dijet
  std::vector<MonitorElement*> dijet_PFJet_pt;
  std::vector<MonitorElement*> dijet_PFJet_eta;
  std::vector<MonitorElement*> dijet_PFJet_phi;
  std::vector<MonitorElement*> dijet_PFJet_rapidity;
  std::vector<MonitorElement*> dijet_PFJet_mass;
  std::vector<MonitorElement*> dijet_deltaPhiPFJet1PFJet2;
  std::vector<MonitorElement*> dijet_deltaEtaPFJet1PFJet2;
  std::vector<MonitorElement*> dijet_deltaRPFJet1PFJet2;
  std::vector<MonitorElement*> dijet_invMassPFJet1PFJet2;
  std::vector<MonitorElement*> dijet_PFchef;
  std::vector<MonitorElement*> dijet_PFnhef;
  std::vector<MonitorElement*> dijet_PFcemf;
  std::vector<MonitorElement*> dijet_PFnemf;
  std::vector<MonitorElement*> dijet_PFJetMulti;
  //
  double dijet_PFJet1_pt_cut_;
  double dijet_PFJet2_pt_cut_;
  int dijet_countPFJet_;

  ///////////////////////////
  // Histograms - DiMuon
  //
  MonitorElement* dimuon_Muon_pt;
  MonitorElement* dimuon_Muon_eta;
  MonitorElement* dimuon_Muon_phi;
  MonitorElement* dimuon_Charge;
  MonitorElement* dimuon_deltaEtaMuon1Muon2;
  MonitorElement* dimuon_deltaPhiMuon1Muon2;
  MonitorElement* dimuon_deltaRMuon1Muon2;
  MonitorElement* dimuon_invMassMuon1Muon2;
  MonitorElement* dimuon_MuonMulti;
  //
  double dimuon_Muon1_pt_cut_;
  double dimuon_Muon2_pt_cut_;
  int dimuon_countMuon_;

  ///////////////////////////
  // Histograms - DiElectron
  //
  MonitorElement* dielectron_Electron_pt;
  MonitorElement* dielectron_Electron_eta;
  MonitorElement* dielectron_Electron_phi;
  MonitorElement* dielectron_Charge;
  MonitorElement* dielectron_deltaEtaElectron1Electron2;
  MonitorElement* dielectron_deltaPhiElectron1Electron2;
  MonitorElement* dielectron_deltaRElectron1Electron2;
  MonitorElement* dielectron_invMassElectron1Electron2;
  MonitorElement* dielectron_ElectronMulti;
  //
  double dielectron_Electron1_pt_cut_;
  double dielectron_Electron2_pt_cut_;
  int dielectron_countElectron_;

  ///////////////////////////
  // Histograms - DiPhoton
  //
  MonitorElement* diphoton_Photon_pt;
  MonitorElement* diphoton_Photon_energy;
  MonitorElement* diphoton_Photon_et;
  MonitorElement* diphoton_Photon_eta;
  MonitorElement* diphoton_Photon_etasc;
  MonitorElement* diphoton_Photon_phi;
  MonitorElement* diphoton_Photon_hovere_eb;
  MonitorElement* diphoton_Photon_hovere_ee;
  MonitorElement* diphoton_Photon_sigmaietaieta_eb;
  MonitorElement* diphoton_Photon_sigmaietaieta_ee;
  MonitorElement* diphoton_Photon_trksumptsolidconedr03_eb;
  MonitorElement* diphoton_Photon_trksumptsolidconedr03_ee;
  MonitorElement* diphoton_Photon_e1x5e5x5_eb;
  MonitorElement* diphoton_Photon_e1x5e5x5_ee;
  MonitorElement* diphoton_Photon_e2x5e5x5_eb;
  MonitorElement* diphoton_Photon_e2x5e5x5_ee;
  MonitorElement* diphoton_deltaEtaPhoton1Photon2;
  MonitorElement* diphoton_deltaPhiPhoton1Photon2;
  MonitorElement* diphoton_deltaRPhoton1Photon2;
  MonitorElement* diphoton_invMassPhoton1Photon2;
  MonitorElement* diphoton_PhotonMulti;
  //
  double diphoton_Photon1_pt_cut_;
  double diphoton_Photon2_pt_cut_;
  int diphoton_countPhoton_;

  ///////////////////////////
  // Histograms - MonoJet
  //
  std::vector<MonitorElement*> monojet_PFJet_pt;
  std::vector<MonitorElement*> monojet_PFJet_eta;
  std::vector<MonitorElement*> monojet_PFJet_phi;
  std::vector<MonitorElement*> monojet_PFMet;
  std::vector<MonitorElement*> monojet_PFMet_phi;
  std::vector<MonitorElement*> monojet_PFJetPtOverPFMet;
  std::vector<MonitorElement*> monojet_deltaPhiPFJetPFMet;
  std::vector<MonitorElement*> monojet_PFchef;
  std::vector<MonitorElement*> monojet_PFnhef;
  std::vector<MonitorElement*> monojet_PFcemf;
  std::vector<MonitorElement*> monojet_PFnemf;
  std::vector<MonitorElement*> monojet_PFJetMulti;
  //
  double monojet_PFJet_pt_cut_;
  double monojet_PFJet_met_cut_;
  int monojet_countPFJet_;

  ///////////////////////////
  // Histograms - MonoMuon
  //
  MonitorElement* monomuon_Muon_pt;
  MonitorElement* monomuon_Muon_eta;
  MonitorElement* monomuon_Muon_phi;
  MonitorElement* monomuon_Charge;
  MonitorElement* monomuon_PFMet;
  MonitorElement* monomuon_PFMet_phi;
  MonitorElement* monomuon_MuonPtOverPFMet;
  MonitorElement* monomuon_deltaPhiMuonPFMet;
  MonitorElement* monomuon_TransverseMass;
  MonitorElement* monomuon_MuonMulti;
  //
  double monomuon_Muon_pt_cut_;
  double monomuon_Muon_met_cut_;
  int monomuon_countMuon_;

  /////////////////////////////
  // Histograms - MonoElectron
  //
  MonitorElement* monoelectron_Electron_pt;
  MonitorElement* monoelectron_Electron_eta;
  MonitorElement* monoelectron_Electron_phi;
  MonitorElement* monoelectron_Charge;
  MonitorElement* monoelectron_PFMet;
  MonitorElement* monoelectron_ElectronPtOverPFMet;
  MonitorElement* monoelectron_PFMet_phi;
  MonitorElement* monoelectron_deltaPhiElectronPFMet;
  MonitorElement* monoelectron_TransverseMass;
  MonitorElement* monoelectron_ElectronMulti;
  //
  double monoelectron_Electron_pt_cut_;
  double monoelectron_Electron_met_cut_;
  int monoelectron_countElectron_;


  ///////////////////////////////////
  // Histograms - Displaced Leptons or Jets
  //
  MonitorElement* dispFerm_isotrack_dxy;
  MonitorElement* dispFerm_electron_dxy;
  MonitorElement* dispFerm_muon_dxy;
  MonitorElement* dispFerm_dsa_dxy;
  MonitorElement* dispFerm_dgl_dxy;
  MonitorElement* dispElec_track_effi_lxy;
  MonitorElement* dispElec_elec_effi_lxy;
  MonitorElement* dispMuon_track_effi_lxy;
  MonitorElement* dispMuon_muon_effi_lxy;
  MonitorElement* dispMuon_muonDisp_effi_lxy;
  MonitorElement* dispMuon_muonDispSA_effi_lxy;
  MonitorElement* dispJet_track_effi_lxy;

  double dispFermion_eta_cut_;
  double dispFermion_pt_cut_;
};

#endif

/* Local Variables: */
/* show-trailing-whitespace: t */
/* truncate-lines: t */
/* End: */
