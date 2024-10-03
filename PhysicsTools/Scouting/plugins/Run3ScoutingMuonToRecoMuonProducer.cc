// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "fastjet/contrib/SoftKiller.hh"

class Run3ScoutingMuonToRecoMuonProducer : public edm::stream::EDProducer<> {
public:
  explicit Run3ScoutingMuonToRecoMuonProducer(const edm::ParameterSet &);
  ~Run3ScoutingMuonToRecoMuonProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
  void endStream() override {}

  /*
  void createPFCandidates(edm::Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle,
                          std::unique_ptr<reco::PFCandidateCollection> &pfcands);
  void createPFCandidatesSK(edm::Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle,
                            std::unique_ptr<reco::PFCandidateCollection> &pfcands);
  reco::PFCandidate createPFCand(Run3ScoutingParticle scoutingparticle);
  void clearVars();
  */
  reco::Track makeRecoTrack(const Run3ScoutingMuon& muon);

private:
  const edm::InputTag srcMuons_;
  const edm::EDGetTokenT<Run3ScoutingMuonCollection> srcMuonToken_;
};

//
// constructors and destructor
//
Run3ScoutingMuonToRecoMuonProducer::Run3ScoutingMuonToRecoMuonProducer(
    const edm::ParameterSet &iConfig)
    : srcMuons_(iConfig.getParameter<edm::InputTag>("srcMuons")),
      srcMuonToken_(consumes<Run3ScoutingMuonCollection>(srcMuons_)) {
  //register products
  produces<reco::MuonCollection>();
  produces<reco::TrackCollection>();
}

Run3ScoutingMuonToRecoMuonProducer::~Run3ScoutingMuonToRecoMuonProducer() = default;

/*
reco::PFCandidate Run3ScoutingMuonToRecoMuonProducer::createPFCand(Run3ScoutingParticle scoutingparticle) {
  auto m = pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
               ? pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->mass()
               : -99.f;
  auto q = pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
               ? pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->charge()
               : -99.f;
  if (m < -90 or q < -90) {
    LogDebug("createPFCand") << "<Run3ScoutingMuonToRecoMuonProducer::createPFCand>:" << std::endl
                             << "Unrecognisable pdgId - skipping particle" << std::endl;
    return reco::PFCandidate();
  }

  float px = scoutingparticle.pt() * cos(scoutingparticle.phi());
  float py = scoutingparticle.pt() * sin(scoutingparticle.phi());
  float pz = scoutingparticle.pt() * sinh(scoutingparticle.eta());
  float p = scoutingparticle.pt() * cosh(scoutingparticle.eta());
  float energy = std::sqrt(p * p + m * m);
  reco::Particle::LorentzVector p4(px, py, pz, energy);

  static const reco::PFCandidate dummy;
  auto pfcand = reco::PFCandidate(q, p4, dummy.translatePdgIdToType(scoutingparticle.pdgId()));

  bool relativeTrackVars = scoutingparticle.relative_trk_vars();
  vertexIndex_.push_back(scoutingparticle.vertex());
  normchi2_.push_back(scoutingparticle.normchi2());
  dz_.push_back(scoutingparticle.dz());
  dxy_.push_back(scoutingparticle.dxy());
  dzsig_.push_back(scoutingparticle.dzsig());
  dxysig_.push_back(scoutingparticle.dxysig());
  lostInnerHits_.push_back(scoutingparticle.lostInnerHits());
  quality_.push_back(scoutingparticle.quality());
  trkPt_.push_back(relativeTrackVars ? scoutingparticle.trk_pt() + scoutingparticle.pt() : scoutingparticle.trk_pt());
  trkEta_.push_back(relativeTrackVars ? scoutingparticle.trk_eta() + scoutingparticle.eta()
                                      : scoutingparticle.trk_eta());
  trkPhi_.push_back(relativeTrackVars ? scoutingparticle.trk_phi() + scoutingparticle.phi()
                                      : scoutingparticle.trk_phi());

  return pfcand;
}
*/

reco::Track Run3ScoutingMuonToRecoMuonProducer::makeRecoTrack(const Run3ScoutingMuon& muon) {

  reco::Track::Vector p3(math::RhoEtaPhiVector(muon.trk_pt(), muon.trk_eta(), muon.trk_phi()));
  reco::Track::Point vtx(muon.trk_vx(), muon.trk_vy(), muon.trk_vz());

  reco::TrackBase::CovarianceMatrix cov;
  cov(0, 0) = pow(muon.trk_qoverpError(), 2);
  cov(0, 1) = muon.trk_qoverp_lambda_cov();
  cov(0, 2) = muon.trk_qoverp_phi_cov();
  cov(0, 3) = muon.trk_qoverp_dxy_cov();
  cov(0, 4) = muon.trk_qoverp_dsz_cov();
  cov(1, 1) = pow(muon.trk_lambdaError(), 2);
  cov(1, 2) = muon.trk_lambda_phi_cov();
  cov(1, 3) = muon.trk_lambda_dxy_cov();
  cov(1, 4) = muon.trk_lambda_dsz_cov();
  cov(2, 2) = pow(muon.trk_phiError(), 2);
  cov(2, 3) = muon.trk_phi_dxy_cov();
  cov(2, 4) = muon.trk_phi_dsz_cov();
  cov(3, 3) = pow(muon.trk_dxyError(), 2);
  cov(3, 4) = muon.trk_dxy_dsz_cov();
  cov(4, 4) = pow(muon.trk_dszError(), 2);

  reco::Track track(muon.trk_chi2(), muon.trk_ndof(), vtx, p3, muon.charge(), cov);

  // Build hit pattern:
  //Run3ScoutingHitPatternPOD hitPattern = muon.trk_hitPattern();
  //reco::HitPattern theHitPattern(hitPattern);
  //for (int i = 0; i < theHitPattern.numberOfAllHits(reco::HitPattern::TRACK_HITS); i++) {
  //  uint16_t hit = theHitPattern.getHitPattern(reco::HitPattern::TRACK_HITS, i);
  //}
  return track;
}

// ------------ method called to produce the data  ------------
void Run3ScoutingMuonToRecoMuonProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup) {
  using namespace edm;

  Handle<Run3ScoutingMuonCollection> srcMuons;
  iEvent.getByToken(srcMuonToken_, srcMuons);

  auto outMuons = std::make_unique<reco::MuonCollection>();
  auto outTracks = std::make_unique<reco::TrackCollection>();

  // Create tracks
  for (unsigned int i = 0; i < srcMuons->size(); i++) {
    const auto& imuon(srcMuons->at(i));
    reco::Track track = makeRecoTrack(imuon);
    outTracks->push_back(track);
  }
  auto outTracksHandle = iEvent.put(std::move(outTracks));

  // Create muons
  for (unsigned int i = 0; i < srcMuons->size(); i++) {
    const auto& imuon(srcMuons->at(i));
    reco::Candidate::PolarLorentzVector p4(imuon.pt(), imuon.eta(), imuon.phi(), imuon.m());
    reco::TrackRef trackRef(outTracksHandle, i);

    reco::Muon muon(trackRef->charge(), reco::Candidate::LorentzVector(p4), trackRef->vertex());
    muon.setType(imuon.type());

    // Only one track per muon in scouting so some arbitration is needed
    if (muon.isGlobalMuon()) {
      muon.setGlobalTrack(trackRef);
      muon.setInnerTrack(trackRef); // Inner track same as L3 track
      muon.setOuterTrack(trackRef); // Outher track same as L3 track
      muon.setBestTrack(reco::Muon::CombinedTrack); // Use L3 track as default
      muon.setTunePBestTrack(reco::Muon::CombinedTrack); // Use L3 track as default
    } else if (muon.isTrackerMuon()) {
      muon.setInnerTrack(trackRef); // Tracker track same as L3 track
      muon.setBestTrack(reco::Muon::InnerTrack); // Use L3 track as default
      muon.setTunePBestTrack(reco::Muon::CombinedTrack); // Use L3 track as default
    }

    outMuons->push_back(muon);
  }

  iEvent.put(std::move(outMuons));

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingMuonToRecoMuonProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcMuons", edm::InputTag("hltScoutingMuonPackerVtx"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingMuonToRecoMuonProducer);
