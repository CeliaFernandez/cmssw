#include "DQM/Physics/src/ExoticaDQM_miniAOD.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace trigger;

typedef vector<string> vstring;

//
// -- Constructor
//
ExoticaDQM_miniAOD::ExoticaDQM_miniAOD(const edm::ParameterSet& ps) {
  edm::LogInfo("ExoticaDQM_miniAOD") << " Starting ExoticaDQM_miniAOD "
                             << "\n";

  typedef std::vector<edm::InputTag> vtag;

  // Get parameters from configuration file
  // Trigger
  TriggerToken_ = consumes<TriggerResults>(ps.getParameter<edm::InputTag>("TriggerResults"));
  HltPaths_  = ps.getParameter<vector<string> >("HltPaths");
  // Vertices
  VertexToken_ = consumes<reco::VertexCollection>(ps.getParameter<InputTag>("vertexCollection"));
  // Electrons
  ElectronToken_ = consumes<std::vector<pat::Electron> >(ps.getParameter<InputTag>("electronCollection"));
  // Muons
  MuonToken_ = consumes<std::vector<pat::Muon> >(ps.getParameter<InputTag>("muonCollection"));
  // Photons
  PhotonToken_ = consumes<std::vector<pat::Photon> >(ps.getParameter<InputTag>("photonCollection"));
  // DiJets
  DiJetPFJetCollection_ = ps.getParameter<std::vector<edm::InputTag> >("DiJetPFJetCollection");
  for (std::vector<edm::InputTag>::const_iterator jetlabel = DiJetPFJetCollection_.begin(),
                                                  jetlabelEnd = DiJetPFJetCollection_.end();
       jetlabel != jetlabelEnd;
       ++jetlabel) {
    DiJetPFJetToken_.push_back(consumes<std::vector<pat::Jet> >(*jetlabel));
  }
  // MonoJet
  MonoJetPFJetCollection_ = ps.getParameter<std::vector<edm::InputTag> >("MonoJetPFJetCollection");
  for (std::vector<edm::InputTag>::const_iterator jetlabel = MonoJetPFJetCollection_.begin(),
                                                  jetlabelEnd = MonoJetPFJetCollection_.end();
       jetlabel != jetlabelEnd;
       ++jetlabel) {
    MonoJetPFJetToken_.push_back(consumes<std::vector<pat::Jet> >(*jetlabel));
  }
  // MET
  PFMETToken_ = consumes<std::vector<pat::MET> >(ps.getParameter<InputTag>("pfMETCollection"));
  // Tracks
  TrackToken_ = consumes<std::vector<pat::IsolatedTrack> >(ps.getParameter<InputTag>("trackCollection"));
  // Displaced collections
  MuonDispSAToken_ = consumes<reco::TrackCollection>(ps.getParameter<InputTag>("displacedSAMuonCollection"));
  // MC Truth
  GenParticleToken_ = consumes<reco::GenParticleCollection>(ps.getParameter<InputTag>("genParticleCollection"));
  //
  magFieldToken_ = esConsumes();

  //Varibles and Cuts for each Module:
  //Dijet
  dijet_PFJet1_pt_cut_ = ps.getParameter<double>("dijet_PFJet1_pt_cut");
  dijet_PFJet2_pt_cut_ = ps.getParameter<double>("dijet_PFJet2_pt_cut");
  //DiMuon
  dimuon_Muon1_pt_cut_ = ps.getParameter<double>("dimuon_Muon1_pt_cut");
  dimuon_Muon2_pt_cut_ = ps.getParameter<double>("dimuon_Muon2_pt_cut");
  //DiElectron
  dielectron_Electron1_pt_cut_ = ps.getParameter<double>("dielectron_Electron1_pt_cut");
  dielectron_Electron2_pt_cut_ = ps.getParameter<double>("dielectron_Electron2_pt_cut");
  //DiPhoton
  diphoton_Photon1_pt_cut_ = ps.getParameter<double>("diphoton_Photon1_pt_cut");
  diphoton_Photon2_pt_cut_ = ps.getParameter<double>("diphoton_Photon2_pt_cut");
  //MonoJet
  monojet_PFJet_pt_cut_ = ps.getParameter<double>("monojet_PFJet_pt_cut");
  monojet_PFJet_met_cut_ = ps.getParameter<double>("monojet_PFJet_met_cut");
  //MonoMuon
  monomuon_Muon_pt_cut_ = ps.getParameter<double>("monomuon_Muon_pt_cut");
  monomuon_Muon_met_cut_ = ps.getParameter<double>("monomuon_Muon_met_cut");
  //MonoElectron
  monoelectron_Electron_pt_cut_ = ps.getParameter<double>("monoelectron_Electron_pt_cut");
  monoelectron_Electron_met_cut_ = ps.getParameter<double>("monoelectron_Electron_met_cut");
  // Displaced lepton or jet
  dispFermion_eta_cut_ = ps.getParameter<double>("dispFermion_eta_cut");
  dispFermion_pt_cut_ = ps.getParameter<double>("dispFermion_pt_cut");
}

//
// -- Destructor
//
ExoticaDQM_miniAOD::~ExoticaDQM_miniAOD() {
  edm::LogInfo("ExoticaDQM_miniAOD") << " Deleting ExoticaDQM "
                             << "\n";
}

//
//  -- Book histograms
//
void ExoticaDQM_miniAOD::bookHistograms(DQMStore::IBooker& bei, edm::Run const&, edm::EventSetup const&) {
  bei.cd();

  //--- DiJet
  for (unsigned int icoll = 0; icoll < DiJetPFJetCollection_.size(); ++icoll) {
    std::stringstream ss;
    ss << "Physics/Exotica_miniAOD/Dijets/" << DiJetPFJetCollection_[icoll].label();
    bei.setCurrentFolder(ss.str());
    dijet_PFJet_pt.push_back(bei.book1D("dijet_PFJet_pt", "PFJet p_{T} (GeV);PFJet p_{T} (GeV); Counts", 50, 30.0, 5000));
    dijet_PFJet_eta.push_back(bei.book1D("dijet_PFJet_eta", "PFJet eta;PFJet #eta; Counts", 50, -2.5, 2.5));
    dijet_PFJet_phi.push_back(bei.book1D("dijet_PFJet_phi", "PFJet #phi;PFJet #phi; Counts", 50, -3.14, 3.14));
    dijet_PFJet_rapidity.push_back(bei.book1D("dijet_PFJet_rapidity", "PFJet rapidity;PFJet rapidity; Counts", 50, -6.0, 6.0));
    dijet_PFJet_mass.push_back(bei.book1D("dijet_PFJet_mass", "PFJet mass (GeV);PFJet mass (GeV); Counts", 50, 0., 500.));
    dijet_deltaPhiPFJet1PFJet2.push_back(
        bei.book1D("dijet_deltaPhiPFJet1PFJet2", "#Delta#phi(Leading PFJet, Sub PFJet);#Delta#phi(Leading PFJet, Sub PFJet); Counts", 40, 0., 3.15));
    dijet_deltaEtaPFJet1PFJet2.push_back(
        bei.book1D("dijet_deltaEtaPFJet1PFJet2", "#Delta#eta(Leading PFJet, Sub PFJet);#Delta#eta(Leading PFJet, Sub PFJet); Counts", 40, -5., 5.));
    dijet_deltaRPFJet1PFJet2.push_back(
        bei.book1D("dijet_deltaRPFJet1PFJet2", "#DeltaR(Leading PFJet, Sub PFJet);#DeltaR(Leading PFJet, Sub PFJet); Counts", 50, 0., 6.));
    dijet_invMassPFJet1PFJet2.push_back(
        bei.book1D("dijet_invMassPFJet1PFJet2", "Leading PFJet, SubLeading PFJet Invariant mass (GeV);Leading PFJet, SubLeading PFJet Invariant mass (GeV); Counts", 50, 0., 8000.));
    dijet_PFchef.push_back(bei.book1D("dijet_PFchef", "Leading PFJet CHEF;Leading PFJet CHEF; Counts", 50, 0.0, 1.0));
    dijet_PFnhef.push_back(bei.book1D("dijet_PFnhef", "Leading PFJet NHEF;Leading PFJet NHEF; Counts", 50, 0.0, 1.0));
    dijet_PFcemf.push_back(bei.book1D("dijet_PFcemf", "Leading PFJet CEMF;Leading PFJet CEMF; Counts", 50, 0.0, 1.0));
    dijet_PFnemf.push_back(bei.book1D("dijet_PFnemf", "Leading PFJEt NEMF;Leading PFJEt NEMF; Counts", 50, 0.0, 1.0));
    dijet_PFJetMulti.push_back(bei.book1D("dijet_PFJetMulti", "No. of PFJets;No. of PFJets; Counts", 10, 0., 10.));
  }
  //--- DiMuon
  bei.setCurrentFolder("Physics/Exotica_miniAOD/DiMuons");
  dimuon_Muon_pt = bei.book1D("dimuon_Muon_pt", "Muon p_{T} (GeV);Muon p_{T} (GeV); Counts", 50, 30.0, 1000);
  dimuon_Muon_eta = bei.book1D("dimuon_Muon_eta", "Muon #eta;Muon #eta; Counts", 50, -2.5, 2.5);
  dimuon_Muon_phi = bei.book1D("dimuon_Muon_phi", "Muon #phi;Muon #phi; Counts", 50, -3.14, 3.14);
  dimuon_Charge = bei.book1D("dimuon_Charge", "Muon charge;Muon charge; Counts", 10, -5., 5.);
  dimuon_deltaEtaMuon1Muon2 =
      bei.book1D("dimuon_deltaEtaMuon1Muon2", "#Delta#eta(Leading Muon, Sub Muon);#Delta#eta(Leading Muon, Sub Muon); Counts", 40, -5., 5.);
  dimuon_deltaPhiMuon1Muon2 =
      bei.book1D("dimuon_deltaPhiMuon1Muon2", "#Delta#phi(Leading Muon, Sub Muon);#Delta#phi(Leading Muon, Sub Muon); Counts", 40, 0., 3.15);
  dimuon_deltaRMuon1Muon2 = bei.book1D("dimuon_deltaRMuon1Muon2", "#DeltaR(Leading Muon, Sub Muon);#DeltaR(Leading Muon, Sub Muon); Counts", 50, 0., 6.);
  dimuon_invMassMuon1Muon2 =
      bei.book1D("dimuon_invMassMuon1Muon2", "Leading Muon, SubLeading Muon Invariant mass (GeV);Leading Muon, SubLeading Muon Invariant mass (GeV); Counts", 50, 70., 2000.);
  dimuon_MuonMulti = bei.book1D("dimuon_MuonMulti", "No. of Muons;No. of Muons; Counts", 10, 0., 10.);
  //--- DiElectrons
  bei.setCurrentFolder("Physics/Exotica_miniAOD/DiElectrons");
  dielectron_Electron_pt = bei.book1D("dielectron_Electron_pt", "Electron p_{T} (GeV);Electron p_{T} (GeV); Counts", 50, 30.0, 1000);
  dielectron_Electron_eta = bei.book1D("dielectron_Electron_eta", "Electron #eta;Electron #eta; Counts", 50, -2.5, 2.5);
  dielectron_Electron_phi = bei.book1D("dielectron_Electron_phi", "Electron #phi;Electron #phi; Counts", 50, -3.14, 3.14);
  dielectron_Charge = bei.book1D("dielectron_Charge", "Electron charge;Electron charge; Counts", 10, -5., 5.);
  dielectron_deltaEtaElectron1Electron2 =
      bei.book1D("dielectron_deltaEtaElectron1Electron2", "#Delta#eta(Leading Electron, Sub Electron);#Delta#eta(Leading Electron, Sub Electron); Counts", 40, -5., 5.);
  dielectron_deltaPhiElectron1Electron2 =
      bei.book1D("dielectron_deltaPhiElectron1Electron2", "#Delta#phi(Leading Electron, Sub Electron);#Delta#phi(Leading Electron, Sub Electron); Counts", 40, 0., 3.15);
  dielectron_deltaRElectron1Electron2 =
      bei.book1D("dielectron_deltaRElectron1Electron2", "#DeltaR(Leading Electron, Sub Electron);#DeltaR(Leading Electron, Sub Electron); Counts", 50, 0., 6.);
  dielectron_invMassElectron1Electron2 = bei.book1D("dielectron_invMassElectron1Electron2",
                                                    "Leading Electron, SubLeading Electron Invariant mass (GeV);Leading Electron, SubLeading Electron Invariant mass (GeV); Counts",
                                                    50,
                                                    70.,
                                                    2000.);
  dielectron_ElectronMulti = bei.book1D("dielectron_ElectronMulti", "No. of Electrons;No. of Electrons; Counts", 10, 0., 10.);
  //--- DiPhotons
  bei.setCurrentFolder("Physics/Exotica_miniAOD/DiPhotons");
  diphoton_Photon_energy = bei.book1D("diphoton_Photon_energy", "Photon energy (GeV);Photon energy (GeV); Counts", 50, 30.0, 300);
  diphoton_Photon_et = bei.book1D("diphoton_Photon_et", "Photon E_{T} (GeV);Photon E_{T} (GeV); Counts", 50, 30.0, 300);
  diphoton_Photon_pt = bei.book1D("diphoton_Photon_pt", "Photon p_{T} (GeV);Photon p_{T} (GeV); Counts", 50, 30.0, 300);
  diphoton_Photon_eta = bei.book1D("diphoton_Photon_eta", "Photon #eta;Photon #eta; Counts", 50, -2.5, 2.5);
  diphoton_Photon_etasc = bei.book1D("diphoton_Photon_etasc", "Photon #eta (SuperCluster);Photon #eta (SuperCluster); Counts", 50, -2.5, 2.5);
  diphoton_Photon_phi = bei.book1D("diphoton_Photon_phi", "Photon #phi;Photon #phi; Counts", 50, -3.14, 3.14);
  diphoton_Photon_hovere_eb = bei.book1D("diphoton_Photon_hovere_eb", "Photon H/E (EB);Photon H/E (EB); Counts", 50, 0., 0.50);
  diphoton_Photon_hovere_ee = bei.book1D("diphoton_Photon_hovere_ee", "Photon H/E (EE);Photon H/E (EE); Counts", 50, 0., 0.50);
  diphoton_Photon_sigmaietaieta_eb =
      bei.book1D("diphoton_Photon_sigmaietaieta_eb", "Photon #sigma_{i #eta i #eta} (EB);Photon #sigma_{i #eta i #eta} (EB); Counts", 50, 0., 0.03);
  diphoton_Photon_sigmaietaieta_ee =
      bei.book1D("diphoton_Photon_sigmaietaieta_ee", "Photon #sigma_{i #eta i #eta} (EE);Photon #sigma_{i #eta i #eta} (EE); Counts", 50, 0., 0.03);
  diphoton_Photon_trksumptsolidconedr03_eb =
      bei.book1D("diphoton_Photon_trksumptsolidconedr03_eb", "Photon TrkSumPtDr03 (EB);Photon TrkSumPtDr03 (EB); Counts", 50, 0., 15.);
  diphoton_Photon_trksumptsolidconedr03_ee =
      bei.book1D("diphoton_Photon_trksumptsolidconedr03_ee", "Photon TrkSumPtDr03 (EE);Photon TrkSumPtDr03 (EE); Counts", 50, 0., 15.);
  diphoton_Photon_e1x5e5x5_eb = bei.book1D("diphoton_Photon_e1x5e5x5_eb", "Photon E_{1x5}/E_{5x5} (EB);Photon E_{1x5}/E_{5x5} (EB); Counts", 50, 0., 1.);
  diphoton_Photon_e1x5e5x5_ee = bei.book1D("diphoton_Photon_e1x5e5x5_ee", "Photon E_{1x5}/E_{5x5} (EE);Photon E_{1x5}/E_{5x5} (EE); Counts", 50, 0., 1.);
  diphoton_Photon_e2x5e5x5_eb = bei.book1D("diphoton_Photon_e2x5e5x5_eb", "Photon E_{2x5}/E_{5x5} (EB);Photon E_{2x5}/E_{5x5} (EB); Counts", 50, 0., 1.);
  diphoton_Photon_e2x5e5x5_ee = bei.book1D("diphoton_Photon_e2x5e5x5_ee", "Photon E_{2x5}/E_{5x5} (EE);Photon E_{2x5}/E_{5x5} (EE); Counts", 50, 0., 1.);
  diphoton_deltaEtaPhoton1Photon2 =
      bei.book1D("diphoton_deltaEtaPhoton1Photon2", "#Delta#eta(Leading Photon, Sub Photon);#Delta#eta(Leading Photon, Sub Photon); Counts", 40, -5., 5.);
  diphoton_deltaPhiPhoton1Photon2 =
      bei.book1D("diphoton_deltaPhiPhoton1Photon2", "#Delta#phi(Leading Photon, Sub Photon);#Delta#phi(Leading Photon, Sub Photon); Counts", 40, 0., 3.15);
  diphoton_deltaRPhoton1Photon2 =
      bei.book1D("diphoton_deltaRPhoton1Photon2", "#DeltaR(Leading Photon, Sub Photon);#DeltaR(Leading Photon, Sub Photon); Counts", 50, 0., 6.);
  diphoton_invMassPhoton1Photon2 = bei.book1D(
      "diphoton_invMassPhoton1Photon2", "Leading Photon, SubLeading Photon Invariant mass (GeV);Leading Photon, SubLeading Photon Invariant mass (GeV); Counts", 50, 500., 4500.);
  diphoton_PhotonMulti = bei.book1D("diphoton_PhotonMulti", "No. of Photons;No. of Photons; Counts", 10, 0., 10.);
  //--- MonoJet
  for (unsigned int icoll = 0; icoll < MonoJetPFJetCollection_.size(); ++icoll) {
    std::stringstream ss;
    ss << "Physics/Exotica_miniAOD/MonoJet/" << MonoJetPFJetCollection_[icoll].label();
    bei.setCurrentFolder(ss.str());
    monojet_PFJet_pt.push_back(bei.book1D("monojet_PFJet_pt", "PFJet p_{T} (GeV);PFJet p_{T} (GeV); Counts", 50, 30.0, 1000));
    monojet_PFJet_eta.push_back(bei.book1D("monojet_PFJet_eta", "PFJet #eta;PFJet #eta; Counts", 50, -2.5, 2.5));
    monojet_PFJet_phi.push_back(bei.book1D("monojet_PFJet_phi", "PFJet #phi;PFJet #phi; Counts", 50, -3.14, 3.14));
    monojet_PFMet.push_back(bei.book1D("monojet_PFMet", "PFMET p_{T} (GeV);PFMET p_{T} (GeV); Counts", 40, 0.0, 1000));
    monojet_PFMet_phi.push_back(bei.book1D("monojet_PFMet_phi", "#phi(PFMET #phi);#phi(PFMET #phi); Counts", 50, -3.14, 3.14));
    monojet_PFJetPtOverPFMet.push_back(bei.book1D("monojet_PFJetPtOverPFMet", "PFJet p_{T} / PFMET;PFJet p_{T} / PFMET; Counts", 40, 0.0, 5.));
    monojet_deltaPhiPFJetPFMet.push_back(bei.book1D("monojet_deltaPhiPFJetPFMet", "#Delta#phi(PFJet, PFMet);#Delta#phi(PFJet, PFMet); Counts", 40, 0., 3.15));
    monojet_PFchef.push_back(bei.book1D("monojet_PFchef", "PFJet CHEF;PFJet CHEF; Counts", 50, 0.0, 1.0));
    monojet_PFnhef.push_back(bei.book1D("monojet_PFnhef", "PFJet NHEF;PFJet NHEF; Counts", 50, 0.0, 1.0));
    monojet_PFcemf.push_back(bei.book1D("monojet_PFcemf", "PFJet CEMF;PFJet CEMF; Counts", 50, 0.0, 1.0));
    monojet_PFnemf.push_back(bei.book1D("monojet_PFnemf", "PFJet NEMF;PFJet NEMF; Counts", 50, 0.0, 1.0));
    monojet_PFJetMulti.push_back(bei.book1D("monojet_PFJetMulti", "No. of PFJets;No. of PFJets; Counts", 10, 0., 10.));
  }
  //--- MonoMuon
  bei.setCurrentFolder("Physics/Exotica_miniAOD/MonoMuon");
  monomuon_Muon_pt = bei.book1D("monomuon_Muon_pt", "Muon p_{T} (GeV);Muon p_{T} (GeV); Counts", 50, 30.0, 2000);
  monomuon_Muon_eta = bei.book1D("monomuon_Muon_eta", "Muon #eta;Muon #eta; Counts", 50, -2.5, 2.5);
  monomuon_Muon_phi = bei.book1D("monomuon_Muon_phi", "Muon #phi;Muon #phi; Counts", 50, -3.14, 3.14);
  monomuon_Charge = bei.book1D("monomuon_Charge", "Muon charge;Muon charge; Counts", 10, -5., 5.);
  monomuon_PFMet = bei.book1D("monomuon_PFMet", "PFMET p_{T} (GeV);PFMET p_{T} (GeV); Counts", 40, 0.0, 2000);
  monomuon_PFMet_phi = bei.book1D("monomuon_PFMet_phi", "PFMET #phi;PFMET #phi; Counts", 50, -3.14, 3.14);
  monomuon_MuonPtOverPFMet = bei.book1D("monomuon_MuonPtOverPFMet", "Muon p_{T} / PFMET;Muon p_{T} / PFMET; Counts", 40, 0.0, 5.);
  monomuon_deltaPhiMuonPFMet = bei.book1D("monomuon_deltaPhiMuonPFMet", "#Delta#phi(Muon, PFMet);#Delta#phi(Muon, PFMet); Counts", 40, 0., 3.15);
  monomuon_TransverseMass = bei.book1D("monomuon_TransverseMass", "Transverse Mass M_{T} (GeV);Transverse Mass M_{T} (GeV); Counts", 40, 200., 3000.);
  monomuon_MuonMulti = bei.book1D("monomuon_MuonMulti", "No. of Muons;No. of Muons; Counts", 10, 0., 10.);
  //--- MonoElectron
  bei.setCurrentFolder("Physics/Exotica_miniAOD/MonoElectron");
  monoelectron_Electron_pt = bei.book1D("monoelectron_Electron_pt", "Electron p_{T} (GeV);Electron p_{T} (GeV); Counts", 50, 30.0, 2000);
  monoelectron_Electron_eta = bei.book1D("monoelectron_Electron_eta", "Electron #eta;Electron #eta; Counts", 50, -2.5, 2.5);
  monoelectron_Electron_phi = bei.book1D("monoelectron_Electron_phi", "Electron #phi;Electron #phi; Counts", 50, -3.14, 3.14);
  monoelectron_Charge = bei.book1D("monoelectron_Charge", "Electron charge;Electron charge; Counts", 10, -5., 5.);
  monoelectron_PFMet = bei.book1D("monoelectron_PFMet", "PFMET p_{T} (GeV);PFMET p_{T} (GeV); Counts", 40, 0.0, 2000);
  monoelectron_PFMet_phi = bei.book1D("monoelectron_PFMet_phi", "PFMET #phi;PFMET #phi; Counts", 50, -3.14, 3.14);
  monoelectron_ElectronPtOverPFMet =
      bei.book1D("monoelectron_ElectronPtOverPFMet", "Electron p_{T} / PFMet;Electron p_{T} / PFMet; Counts", 40, 0.0, 5.);
  monoelectron_deltaPhiElectronPFMet =
      bei.book1D("monoelectron_deltaPhiElectronPFMet", "#Delta#phi(Electron, PFMet);#Delta#phi(Electron, PFMet); Counts", 40, 0., 3.15);
  monoelectron_TransverseMass = bei.book1D("monoelectron_TransverseMass", "Transverse Mass M_{T} (GeV);Transverse Mass M_{T} (GeV); Counts", 40, 200., 4000.);
  monoelectron_ElectronMulti = bei.book1D("monoelectron_ElectronMulti", "No. of Electrons;No. of Electrons; Counts", 10, 0., 10.);


  //--- Displaced Leptons (filled using only leptons from long-lived stop decay).
  bei.setCurrentFolder("Physics/Exotica_miniAOD/DisplacedFermions");
  dispElec_track_effi_lxy = bei.bookProfile("dispElec_track_effi_lxy",
                                            "Electron channel; Transverse decay length (cm); Track reco efficiency",
                                            10,
                                            0.,
                                            60.,
                                            -999.,
                                            999,
                                            "");
  dispElec_elec_effi_lxy = bei.bookProfile("dispElec_elec_effi_lxy",
                                           "Electron channel; Transverse decay length (cm); Electron reco efficiency",
                                           10,
                                           0.,
                                           60.,
                                           -999.,
                                           999,
                                           "");
  dispMuon_track_effi_lxy = bei.bookProfile("dispMuon_track_effi_lxy",
                                            "Muon channel; Transverse decay length (cm); Track reco efficiency",
                                            10,
                                            0.,
                                            60.,
                                            -999.,
                                            999,
                                            "");
  dispMuon_muon_effi_lxy = bei.bookProfile("dispMuon_muon_effi_lxy",
                                           "Muon channel; Transverse decay length (cm); Muon reco efficiency",
                                           10,
                                           0.,
                                           100.,
                                           -999.,
                                           999,
                                           "");
  dispMuon_muonDisp_effi_lxy =
      bei.bookProfile("dispMuon_muonDisp_effi_lxy",
                      "Muon channel; Transverse decay length (cm); dGlobalMuon reco efficiency",
                      10,
                      0.,
                      100.,
                      -999.,
                      999,
                      "");
  dispMuon_muonDispSA_effi_lxy =
      bei.bookProfile("dispMuon_muonDispSA_effi_lxy",
                      "Muon channel; Transverse decay length (cm); dSAMuon reco efficiency",
                      10,
                      0.,
                      400.,
                      -999.,
                      999,
                      "");
  //--- Displaced Jets (filled using only tracks or jets from long-lived stop decay).
  dispJet_track_effi_lxy = bei.bookProfile("dispJet_track_effi_lxy",
                                           "Jet channel; Transverse decay length (cm); Track reco efficiency",
                                           10,
                                           0.,
                                           100.,
                                           -999.,
                                           999,
                                           "");

  bei.cd();
}

//
//  -- Analyze
//
void ExoticaDQM_miniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // objects

  //Trigger
  bool ValidTriggers = iEvent.getByToken(TriggerToken_, TriggerResults_);
  if (!ValidTriggers) {
    return;
  }

  // Vertices
  bool ValidVertices = iEvent.getByToken(VertexToken_, VertexCollection_);
  if (!ValidVertices) {
    return;
  }

  // Electrons
  iEvent.getByToken(ElectronToken_, ElectronCollection_);
  bool ValidElectrons = iEvent.getByToken(ElectronToken_, ElectronCollection_);
  if (!ValidElectrons) {
    return;
  }

  // Muons
  bool ValidPFMuon = iEvent.getByToken(MuonToken_, MuonCollection_);
  if (!ValidPFMuon) {
    return;
  }

  // PFMETs
  bool ValidPFMET = iEvent.getByToken(PFMETToken_, pfMETCollection_);
  if (!ValidPFMET) {
    return;
  }

  // Photons
  bool ValidCaloPhoton = iEvent.getByToken(PhotonToken_, PhotonCollection_);
  if (!ValidCaloPhoton) {
    return;
  }

  // Tracks
  bool ValidTracks = iEvent.getByToken(TrackToken_, TrackCollection_);

  // Special collections for displaced particles
  bool ValidDSA = iEvent.getByToken(MuonDispSAToken_, MuonDispSACollection_);

  // MC truth
  bool ValidGenParticles = iEvent.getByToken(GenParticleToken_, GenCollection_);

  // ---> Trigger

  int N_Triggers = TriggerResults_->size();
  int N_GoodTriggerPaths = HltPaths_.size();
  bool triggered_event = false;
  const edm::TriggerNames& trigName = iEvent.triggerNames(*TriggerResults_);
  for (int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig) {
    if (TriggerResults_.product()->accept(i_Trig)) {
      for (int n = 0; n < N_GoodTriggerPaths; n++) {
        if (trigName.triggerName(i_Trig).find(HltPaths_[n]) != std::string::npos) {
          triggered_event = true;
        }
      }
    }
  }
  if (triggered_event == false)
    return;

  // ---> Variable initialisation
  for (int i = 0; i < 2; i++) {
    //Jets
    PFJetPx[i] = 0.;
    PFJetPy[i] = 0.;
    PFJetPt[i] = 0.;
    PFJetEta[i] = 0.;
    PFJetPhi[i] = 0.;
    PFJetMass[i] = 0.;
    PFJetNHEF[i] = 0.;
    PFJetCHEF[i] = 0.;
    PFJetNEMF[i] = 0.;
    PFJetCEMF[i] = 0.;
    //Muons
    MuonPx[i] = 0.;
    MuonPy[i] = 0.;
    MuonPt[i] = 0.;
    MuonEta[i] = 0.;
    MuonPhi[i] = 0.;
    MuonCharge[i] = 0.;
    //Electrons
    ElectronPx[i] = 0.;
    ElectronPy[i] = 0.;
    ElectronPt[i] = 0.;
    ElectronEta[i] = 0.;
    ElectronPhi[i] = 0.;
    ElectronCharge[i] = 0.;
    //Photons
    PhotonEnergy[i] = 0.;
    PhotonPt[i] = 0.;
    PhotonEt[i] = 0.;
    PhotonEta[i] = 0.;
    PhotonEtaSc[i] = 0.;
    PhotonPhi[i] = 0.;
    PhotonHoverE[i] = 0.;
    PhotonSigmaIetaIeta[i] = 0.;
    PhotonTrkSumPtSolidConeDR03[i] = 0.;
    PhotonE1x5E5x5[i] = 0.;
    PhotonE2x5E5x5[i] = 0.;
  }

  // ---> Vertex Collection
  VertexCollection vertexCollection = *(VertexCollection_.product());
  reco::VertexCollection::const_iterator primaryVertex_ = vertexCollection.begin();

  // ---> Muon Collection
  dimuon_countMuon_ = 0;
  monomuon_countMuon_ = 0;
  std::vector<pat::Muon>::const_iterator muon_ = MuonCollection_->begin();
  for (; muon_ != MuonCollection_->end(); muon_++) {
    // Muon High Pt ID
    if (muon_->isHighPtMuon(*primaryVertex_) == true) {
      if (muon_->pt() > MuonPt[0]) {
        MuonPt[1] = MuonPt[0];
        MuonPx[1] = MuonPx[0];
        MuonPy[1] = MuonPy[0];
        MuonEta[1] = MuonEta[0];
        MuonPhi[1] = MuonPhi[0];
        MuonCharge[1] = MuonCharge[0];
        //
        MuonPt[0] = muon_->pt();
        MuonPx[0] = muon_->px();
        MuonPy[0] = muon_->py();
        MuonEta[0] = muon_->eta();
        MuonPhi[0] = muon_->phi();
        MuonCharge[0] = muon_->charge();
      } else if (muon_->pt() < MuonPt[0] && muon_->pt() > MuonPt[1]) {
        MuonPt[1] = muon_->pt();
        MuonPx[1] = muon_->px();
        MuonPy[1] = muon_->py();
        MuonEta[1] = muon_->eta();
        MuonPhi[1] = muon_->phi();
        MuonCharge[1] = muon_->charge();
      } else {
      }
    }
    if (muon_->pt() > dimuon_Muon1_pt_cut_)
      dimuon_countMuon_++;
    if (muon_->pt() > dimuon_Muon1_pt_cut_)
      monomuon_countMuon_++;
  }

  // -> Electron Collection
  dielectron_countElectron_ = 0;
  monoelectron_countElectron_ = 0;
  std::vector<pat::Electron>::const_iterator electron_ = ElectronCollection_->begin();
  for (; electron_ != ElectronCollection_->end(); electron_++) {
    //HEEP Selection 4.1 (some cuts)
    if (electron_->e5x5() <= 0)
      continue;
    if (electron_->gsfTrack().isNull())
      continue;
    bool HEPP_ele = false;
    double sceta = electron_->caloPosition().eta();
    double dEtaIn = fabs(electron_->deltaEtaSuperClusterTrackAtVtx());
    double dPhiIn = fabs(electron_->deltaPhiSuperClusterTrackAtVtx());
    double HoverE = electron_->hadronicOverEm();
    int missingHits = electron_->gsfTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
    double dxy = electron_->gsfTrack()->dxy(primaryVertex_->position());
    double tkIso = electron_->dr03TkSumPt();
    double e2x5Fraction = electron_->e2x5Max() / electron_->e5x5();
    double e1x5Fraction = electron_->e1x5() / electron_->e5x5();
    double scSigmaIetaIeta = electron_->scSigmaIEtaIEta();
    if (electron_->ecalDriven() && electron_->pt() > dielectron_Electron1_pt_cut_) {
      if (fabs(sceta) < 1.442) {  // barrel
        if (fabs(dEtaIn) < 0.005 && fabs(dPhiIn) < 0.06 && HoverE < 0.05 && tkIso < 5. && missingHits <= 1 &&
            fabs(dxy) < 0.02 && (e2x5Fraction > 0.94 || e1x5Fraction > 0.83))
          HEPP_ele = true;
      } else if (fabs(sceta) > 1.56 && fabs(sceta) < 2.5) {  // endcap
        if (fabs(dEtaIn) < 0.007 && fabs(dPhiIn) < 0.06 && HoverE < 0.05 && tkIso < 5. && missingHits <= 1 &&
            fabs(dxy) < 0.02 && scSigmaIetaIeta < 0.03)
          HEPP_ele = true;
      }
    }
    //
    if (HEPP_ele == false)
      continue;
    if (electron_->pt() > ElectronPt[0]) {
      ElectronPt[1] = ElectronPt[0];
      ElectronPx[1] = ElectronPx[0];
      ElectronPy[1] = ElectronPy[0];
      ElectronEta[1] = ElectronEta[0];
      ElectronPhi[1] = ElectronPhi[0];
      ElectronCharge[1] = ElectronCharge[0];
      //
      ElectronPt[0] = electron_->pt();
      ElectronPx[0] = electron_->px();
      ElectronPy[0] = electron_->py();
      ElectronEta[0] = electron_->eta();
      ElectronPhi[0] = electron_->phi();
      ElectronCharge[0] = electron_->charge();
    } else if (electron_->pt() < ElectronPt[0] && electron_->pt() > ElectronPt[1]) {
      ElectronPt[1] = electron_->pt();
      ElectronPx[1] = electron_->px();
      ElectronPy[1] = electron_->py();
      ElectronEta[1] = electron_->eta();
      ElectronPhi[1] = electron_->phi();
      ElectronCharge[1] = electron_->charge();
    } else {
    }
    if (electron_->pt() > dielectron_Electron1_pt_cut_)
      dielectron_countElectron_++;
    if (electron_->pt() > dielectron_Electron1_pt_cut_)
      monoelectron_countElectron_++;
  }

  // -> Photon Collection
  diphoton_countPhoton_ = 0.;
  for (std::vector<pat::Photon>::const_iterator photon_ = PhotonCollection_->begin(); photon_ != PhotonCollection_->end(); ++photon_) {
    if (photon_->pt() > PhotonPt[0]) {
      PhotonEnergy[1] = PhotonEnergy[0];
      PhotonPt[1] = PhotonPt[0];
      PhotonEt[1] = PhotonEt[0];
      PhotonEta[1] = PhotonEta[0];
      PhotonEtaSc[1] = PhotonEtaSc[0];
      PhotonPhi[1] = PhotonPhi[0];
      PhotonHoverE[1] = PhotonHoverE[0];
      PhotonSigmaIetaIeta[1] = PhotonSigmaIetaIeta[0];
      PhotonTrkSumPtSolidConeDR03[1] = PhotonTrkSumPtSolidConeDR03[0];
      PhotonE1x5E5x5[1] = PhotonE1x5E5x5[0];
      PhotonE2x5E5x5[1] = PhotonE2x5E5x5[0];
      //
      PhotonEnergy[0] = photon_->energy();
      PhotonPt[0] = photon_->pt();
      PhotonEt[0] = photon_->et();
      PhotonEta[0] = photon_->eta();
      PhotonEtaSc[0] = photon_->caloPosition().eta();
      PhotonPhi[0] = photon_->phi();
      PhotonHoverE[0] = photon_->hadronicOverEm();
      PhotonSigmaIetaIeta[0] = photon_->sigmaIetaIeta();
      PhotonTrkSumPtSolidConeDR03[0] = photon_->trkSumPtSolidConeDR03();
      PhotonE1x5E5x5[0] = photon_->e1x5() / photon_->e5x5();
      PhotonE2x5E5x5[0] = photon_->e2x5() / photon_->e5x5();
    } else if (photon_->pt() < PhotonPt[0] && photon_->pt() > PhotonPt[1]) {
      PhotonEnergy[1] = photon_->energy();
      PhotonPt[1] = photon_->pt();
      PhotonEt[1] = photon_->et();
      PhotonEta[1] = photon_->eta();
      PhotonEtaSc[1] = photon_->caloPosition().eta();
      PhotonPhi[1] = photon_->phi();
      PhotonHoverE[1] = photon_->hadronicOverEm();
      PhotonSigmaIetaIeta[1] = photon_->sigmaIetaIeta();
      PhotonTrkSumPtSolidConeDR03[1] = photon_->trkSumPtSolidConeDR03();
      PhotonE1x5E5x5[1] = photon_->e1x5() / photon_->e5x5();
      PhotonE2x5E5x5[1] = photon_->e2x5() / photon_->e5x5();
    } else {
    }
    if (photon_->pt() > dielectron_Electron1_pt_cut_)
      diphoton_countPhoton_++;
  }

  // ---> MET Collection
  const pat::MET &pfmet = (*pfMETCollection_)[0];
  PFMETEt = pfmet.et();
  PFMETPhi = pfmet.phi();

  //
  // Analyze
  //

  //Resonances
  analyzeDiJets(iEvent);
  analyzeDiMuons(iEvent);
  analyzeDiElectrons(iEvent);
  analyzeDiPhotons(iEvent);

  //MonoSearches
  analyzeMonoJets(iEvent);
  analyzeMonoMuons(iEvent);
  analyzeMonoElectrons(iEvent);

  //LongLived
  if (ValidGenParticles && ValidTracks && ValidDSA) {
    analyzeDisplacedLeptons(iEvent, iSetup);
    analyzeDisplacedJets(iEvent, iSetup);
  }


}

void ExoticaDQM_miniAOD::analyzeDisplacedLeptons(const Event& iEvent, const edm::EventSetup& iSetup) {
  //=== This is designed to run on MC events in which a pair of long-lived stop quarks each decay to a displaced lepton + displaced b jet.

  // Initialisation

  const unsigned int stop1 = 1000006;  // PDG identifier of top squark1
  const unsigned int stop2 = 2000006;  // PDG identifier of top squark2
  const float deltaRcut = 0.01;        // Cone size for matching reco to true leptons.
  const float invPtcut = 0.1;          // Cut in 1/Pt consistency for matching reco tracks to genParticles.

  //--- Measure the efficiency to reconstruct leptons from long-lived stop quark decay.

  for (const reco::GenParticle& gen : *GenCollection_) {
    unsigned int idPdg = abs(gen.pdgId());
    // Find electrons/muons from long-lived stop decay.
    if (idPdg == stop1 || idPdg == stop2) {
      unsigned int nDau = gen.numberOfDaughters();
      for (unsigned int i = 0; i < nDau; i++) {
        const reco::GenParticle* dau = (const reco::GenParticle*)gen.daughter(i);
        // Only measure efficiency using leptons passing pt & eta cuts. (The pt cut is almost irrelevant, since leptons from stop decay are hard).
        if (fabs(dau->eta()) < dispFermion_eta_cut_ && dau->pt() > dispFermion_pt_cut_) {
          unsigned int pdgIdDau = abs(dau->pdgId());

          if (pdgIdDau == 11 || pdgIdDau == 13) {  // electron or muon from stop decay

            // Get transverse decay length of stop quark.
            float lxy = dau->vertex().rho();

            // Get momentum vector of daughter genParticle trajectory extrapolated to beam-line.
            GlobalVector genP = this->getGenParticleTrajectoryAtBeamline(iSetup, dau);

            if (pdgIdDau == 11) {  // electron from stop decay

              // Find matching reco track if any.
              bool recoedTrk = false;
              for (const pat::IsolatedTrack& trk : *TrackCollection_) {
                if (reco::deltaR(genP, trk) < deltaRcut && fabs(1 / dau->pt() - 1 / trk.pt()) < invPtcut) {
                  recoedTrk = true;
                }
              }
              dispElec_track_effi_lxy->Fill(lxy, recoedTrk);

              // Find matching reco electron if any.
              bool recoedE = false;
              for (const pat::Electron& eReco : *ElectronCollection_) {
                if (reco::deltaR(genP, eReco) < deltaRcut && fabs(1 / dau->pt() - 1 / eReco.pt()) < invPtcut)
                  recoedE = true;
              }
              dispElec_elec_effi_lxy->Fill(lxy, recoedE);

            } else if (pdgIdDau == 13) {  // muon from stop decay

              // Find matching reco track if any.
              bool recoedTrk = false;
              for (const pat::IsolatedTrack& trk : *TrackCollection_) {
                if (reco::deltaR(genP, trk) < deltaRcut && fabs(1 / dau->pt() - 1 / trk.pt()) < invPtcut) {
                    recoedTrk = true;
                }
              }
              dispMuon_track_effi_lxy->Fill(lxy, recoedTrk);

              // Find matching reco muon, if any, in normal muon collection.
              bool recoedMu = false;
              for (const pat::Muon& muReco : *MuonCollection_) {
                if (reco::deltaR(genP, muReco) < deltaRcut && fabs(1 / dau->pt() - 1 / muReco.pt()) < invPtcut)
                  recoedMu = true;
              }
              dispMuon_muon_effi_lxy->Fill(lxy, recoedMu);

              // Find matching reco muon, if any, in displaced SA muon collection.
              bool recoedMuDispSA = false;
              for (const reco::Track& muDispSAReco : *MuonDispSACollection_) {
                if (reco::deltaR(genP, muDispSAReco) < deltaRcut &&
                    fabs(1 / dau->pt() - 1 / muDispSAReco.pt()) < invPtcut)
                  recoedMuDispSA = true;
              }
              dispMuon_muonDispSA_effi_lxy->Fill(lxy, recoedMuDispSA);
            }
          }
        }
      }
    }
  }
}
void ExoticaDQM_miniAOD::analyzeDisplacedJets(const Event& iEvent, const edm::EventSetup& iSetup) {
  //=== This is designed to run on MC events in which a pair of long-lived stop quarks each decay to a displaced lepton + displaced b jet.

  // Initialisation

  // Define function to identify R-hadrons containing stop quarks from PDG particle code.
  // N.B. Jets originate not just from stop quark, but also from its partner SM quark inside the R hadron.
  auto isRhadron = [](unsigned int pdgId) { return (pdgId / 100) == 10006 || (pdgId / 1000) == 1006; };

  const float deltaRcut = 0.01;  // Cone size for matching reco tracks to genParticles.
  const float invPtcut = 0.1;    // Cut in 1/Pt consistency for matching reco tracks to genParticles.

  //--- Measure the efficiency to reconstruct tracks in jet(s) from long-lived stop quark decay.

  for (const reco::GenParticle& gen : *GenCollection_) {
    unsigned int idPdg = abs(gen.pdgId());
    // Only measure efficiency using charged e, mu pi, K, p
    if (idPdg == 11 || idPdg == 13 || idPdg == 211 || idPdg == 321 || idPdg == 2212) {
      // Only measure efficiency using leptons passing pt & eta cuts. (The pt cut is almost irrelevant, since leptons from stop decay are hard).
      if (fabs(gen.eta()) < dispFermion_eta_cut_ && gen.pt() > dispFermion_pt_cut_) {
        // Check if this particle came (maybe indirectly) from an R hadron decay.
        const reco::GenParticle* genMoth = &gen;
        const reco::GenParticle* genRhadron = nullptr;
        bool foundParton = false;
        while (genMoth->numberOfMothers() > 0) {
          genMoth = (const reco::GenParticle*)genMoth->mother(0);
          unsigned int idPdgMoth = abs(genMoth->pdgId());
          // Check that the R-hadron decayed via a quark/gluon before yielding genParticle "gen".
          // This ensures that gen is from the jet, and not a lepton produced directly from the stop quark decay.
          if ((idPdgMoth >= 1 && idPdgMoth <= 6) || idPdgMoth == 21)
            foundParton = true;
          // Note if ancestor was R hadron
          if (isRhadron(idPdgMoth)) {
            genRhadron = genMoth;
            break;
          }
        }

        if (foundParton && genRhadron != nullptr) {  // This GenParticle came (maybe indirectly) from an R hadron decay.

          // Get transverse decay length of R hadron.
          float lxy = genRhadron->daughter(0)->vertex().rho();

          // Get momentum vector of genParticle trajectory extrapolated to beam-line.
          GlobalVector genP = this->getGenParticleTrajectoryAtBeamline(iSetup, &gen);

          // Find matching reco track if any.
          bool recoedTrk = false;
          for (const pat::IsolatedTrack& trk : *TrackCollection_) {
            if (reco::deltaR(genP, trk) < deltaRcut && fabs(1 / gen.pt() - 1 / trk.pt()) < invPtcut) {
              //cout<<"MATCH TRK "<<gen.pt()<<" "<<trk.pt()<<" "<<reco::deltaR(gen, trk)<<endl;
              recoedTrk = true;
            }
          }
          dispJet_track_effi_lxy->Fill(lxy, recoedTrk);
        }
      }
    }
  }
}
GlobalVector ExoticaDQM_miniAOD::getGenParticleTrajectoryAtBeamline(const edm::EventSetup& iSetup,
                                                            const reco::GenParticle* gen) {
  //=== Estimate the momentum vector that a GenParticle would have at its trajectory's point of closest
  //=== approach to the beam-line.

  // Get the magnetic field
  const MagneticField* theMagField = &iSetup.getData(magFieldToken_);

  // Make FreeTrajectoryState of this gen particle
  FreeTrajectoryState fts(GlobalPoint(gen->vx(), gen->vy(), gen->vz()),
                          GlobalVector(gen->px(), gen->py(), gen->pz()),
                          gen->charge(),
                          theMagField);

  // Get trajectory closest to beam line
  TSCBLBuilderNoMaterial tscblBuilder;
  const BeamSpot beamspot;  // Simple beam-spot at (0,0,0). Good enough.
  TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(fts, beamspot);

  GlobalVector p = tsAtClosestApproach.trackStateAtPCA().momentum();

  return p;
}

void ExoticaDQM_miniAOD::analyzeDiJets(const Event& iEvent) {
  for (unsigned int icoll = 0; icoll < DiJetPFJetCollection_.size(); ++icoll) {
    dijet_countPFJet_ = 0;
    bool ValidDiJetPFJets = iEvent.getByToken(DiJetPFJetToken_[icoll], DiJetpfJetCollection_);
    if (!ValidDiJetPFJets)
      continue;
    DiJetpfjets = *DiJetpfJetCollection_;
    for (int i = 0; i < 2; i++) {
      PFJetPx[i] = 0.;
      PFJetPy[i] = 0.;
      PFJetPt[i] = 0.;
      PFJetEta[i] = 0.;
      PFJetPhi[i] = 0.;
      PFJetNHEF[i] = 0.;
      PFJetCHEF[i] = 0.;
      PFJetNEMF[i] = 0.;
      PFJetCEMF[i] = 0.;
    }
    std::vector<pat::Jet>::const_iterator DiJetpfjet_ = DiJetpfjets.begin();
    for (; DiJetpfjet_ != DiJetpfjets.end(); ++DiJetpfjet_) {
      if (!DiJetpfjet_->isPFJet()) { continue; }
      if (DiJetpfjet_->pt() > PFJetPt[0]) {
        PFJetPt[1] = PFJetPt[0];
        PFJetPx[1] = PFJetPx[0];
        PFJetPy[1] = PFJetPy[0];
        PFJetEta[1] = PFJetEta[0];
        PFJetPhi[1] = PFJetPhi[0];
        PFJetRapidity[1] = DiJetpfjet_->rapidity();
        PFJetMass[1] = DiJetpfjet_->mass();
        PFJetNHEF[1] = PFJetNHEF[0];
        PFJetCHEF[1] = PFJetCHEF[0];
        PFJetNEMF[1] = PFJetNEMF[0];
        PFJetCEMF[1] = PFJetCEMF[0];
        PFJetPt[0] = DiJetpfjet_->pt();
        PFJetPx[0] = DiJetpfjet_->px();
        PFJetPy[0] = DiJetpfjet_->py();
        PFJetEta[0] = DiJetpfjet_->eta();
        PFJetPhi[0] = DiJetpfjet_->phi();
        PFJetRapidity[0] = DiJetpfjet_->rapidity();
        PFJetMass[0] = DiJetpfjet_->mass();
        PFJetNHEF[0] = DiJetpfjet_->neutralHadronEnergyFraction();
        PFJetCHEF[0] = DiJetpfjet_->chargedHadronEnergyFraction();
        PFJetNEMF[0] = DiJetpfjet_->neutralEmEnergyFraction();
        PFJetCEMF[0] = DiJetpfjet_->chargedEmEnergyFraction();
      } else if (DiJetpfjet_->pt() < PFJetPt[0] && DiJetpfjet_->pt() > PFJetPt[1]) {
        PFJetPt[1] = DiJetpfjet_->pt();
        PFJetPx[1] = DiJetpfjet_->px();
        PFJetPy[1] = DiJetpfjet_->py();
        PFJetEta[1] = DiJetpfjet_->eta();
        PFJetPhi[1] = DiJetpfjet_->phi();
        PFJetRapidity[1] = DiJetpfjet_->rapidity();
        PFJetMass[1] = DiJetpfjet_->mass();
        PFJetNHEF[1] = DiJetpfjet_->neutralHadronEnergyFraction();
        PFJetCHEF[1] = DiJetpfjet_->chargedHadronEnergyFraction();
        PFJetNEMF[1] = DiJetpfjet_->neutralEmEnergyFraction();
        PFJetCEMF[1] = DiJetpfjet_->chargedEmEnergyFraction();
      } else {
      }
      if (DiJetpfjet_->pt() > dijet_PFJet1_pt_cut_)
        dijet_countPFJet_++;
    }
    if (PFJetPt[0] > dijet_PFJet1_pt_cut_ && PFJetPt[1] > dijet_PFJet2_pt_cut_) {
      dijet_PFJet_pt[icoll]->Fill(PFJetPt[0]);
      dijet_PFJet_eta[icoll]->Fill(PFJetEta[0]);
      dijet_PFJet_phi[icoll]->Fill(PFJetPhi[0]);
      dijet_PFJet_rapidity[icoll]->Fill(PFJetRapidity[0]);
      dijet_PFJet_mass[icoll]->Fill(PFJetMass[0]);
      dijet_PFJet_pt[icoll]->Fill(PFJetPt[1]);
      dijet_PFJet_eta[icoll]->Fill(PFJetEta[1]);
      dijet_PFJet_phi[icoll]->Fill(PFJetPhi[1]);
      dijet_PFJet_rapidity[icoll]->Fill(PFJetRapidity[1]);
      dijet_PFJet_mass[icoll]->Fill(PFJetMass[1]);
      dijet_deltaPhiPFJet1PFJet2[icoll]->Fill(deltaPhi(PFJetPhi[0], PFJetPhi[1]));
      dijet_deltaEtaPFJet1PFJet2[icoll]->Fill(PFJetEta[0] - PFJetEta[1]);
      dijet_deltaRPFJet1PFJet2[icoll]->Fill(deltaR(PFJetEta[0], PFJetPhi[0], PFJetEta[1], PFJetPhi[1]));
      dijet_invMassPFJet1PFJet2[icoll]->Fill(
          sqrt(2 * PFJetPt[0] * PFJetPt[1] * (cosh(PFJetEta[0] - PFJetEta[1]) - cos(PFJetPhi[0] - PFJetPhi[1]))));
      dijet_PFchef[icoll]->Fill(PFJetCHEF[0]);
      dijet_PFnhef[icoll]->Fill(PFJetNHEF[0]);
      dijet_PFcemf[icoll]->Fill(PFJetCEMF[0]);
      dijet_PFnemf[icoll]->Fill(PFJetNEMF[0]);
      dijet_PFchef[icoll]->Fill(PFJetCHEF[1]);
      dijet_PFnhef[icoll]->Fill(PFJetNHEF[1]);
      dijet_PFcemf[icoll]->Fill(PFJetCEMF[1]);
      dijet_PFnemf[icoll]->Fill(PFJetNEMF[1]);
      dijet_PFJetMulti[icoll]->Fill(dijet_countPFJet_);
    }
  }
}

void ExoticaDQM_miniAOD::analyzeDiMuons(const Event& iEvent) {
  if (MuonPt[0] > dimuon_Muon1_pt_cut_ && MuonPt[1] > dimuon_Muon2_pt_cut_ && MuonCharge[0] * MuonCharge[1] == -1) {
    dimuon_Muon_pt->Fill(MuonPt[0]);
    dimuon_Muon_eta->Fill(MuonEta[0]);
    dimuon_Muon_phi->Fill(MuonPhi[0]);
    dimuon_Muon_pt->Fill(MuonPt[1]);
    dimuon_Muon_eta->Fill(MuonEta[1]);
    dimuon_Muon_phi->Fill(MuonPhi[1]);
    dimuon_Charge->Fill(MuonCharge[0]);
    dimuon_Charge->Fill(MuonCharge[1]);
    dimuon_deltaPhiMuon1Muon2->Fill(deltaPhi(MuonPhi[0], MuonPhi[1]));
    dimuon_deltaEtaMuon1Muon2->Fill(MuonEta[0] - MuonEta[1]);
    dimuon_deltaRMuon1Muon2->Fill(deltaR(MuonEta[0], MuonPhi[0], MuonEta[1], MuonPhi[1]));
    dimuon_invMassMuon1Muon2->Fill(
        sqrt(2 * MuonPt[0] * MuonPt[1] * (cosh(MuonEta[0] - MuonEta[1]) - cos(MuonPhi[0] - MuonPhi[1]))));
    dimuon_MuonMulti->Fill(dimuon_countMuon_);
  }
}

void ExoticaDQM_miniAOD::analyzeDiElectrons(const Event& iEvent) {
  if (ElectronPt[0] > dielectron_Electron1_pt_cut_ && ElectronPt[1] > dielectron_Electron2_pt_cut_ &&
      ElectronCharge[0] * ElectronCharge[1] == -1.) {
    dielectron_Electron_pt->Fill(ElectronPt[0]);
    dielectron_Electron_eta->Fill(ElectronEta[0]);
    dielectron_Electron_phi->Fill(ElectronPhi[0]);
    dielectron_Electron_pt->Fill(ElectronPt[1]);
    dielectron_Electron_eta->Fill(ElectronEta[1]);
    dielectron_Electron_phi->Fill(ElectronPhi[1]);
    dielectron_Charge->Fill(ElectronCharge[0]);
    dielectron_Charge->Fill(ElectronCharge[1]);
    dielectron_deltaPhiElectron1Electron2->Fill(deltaPhi(ElectronPhi[0], ElectronPhi[1]));
    dielectron_deltaEtaElectron1Electron2->Fill(ElectronEta[0] - ElectronEta[1]);
    dielectron_deltaRElectron1Electron2->Fill(deltaR(ElectronEta[0], ElectronPhi[0], ElectronEta[1], ElectronPhi[1]));
    dielectron_invMassElectron1Electron2->Fill(
        sqrt(2 * ElectronPt[0] * ElectronPt[1] *
             (cosh(ElectronEta[0] - ElectronEta[1]) - cos(ElectronPhi[0] - ElectronPhi[1]))));
    dielectron_ElectronMulti->Fill(dielectron_countElectron_);
  }
}

void ExoticaDQM_miniAOD::analyzeDiPhotons(const Event& iEvent) {
  if (PhotonPt[0] > diphoton_Photon1_pt_cut_ && PhotonPt[1] > diphoton_Photon2_pt_cut_) {
    diphoton_Photon_energy->Fill(PhotonEnergy[0]);
    diphoton_Photon_pt->Fill(PhotonPt[0]);
    diphoton_Photon_et->Fill(PhotonEt[0]);
    diphoton_Photon_eta->Fill(PhotonEta[0]);
    diphoton_Photon_etasc->Fill(PhotonEtaSc[0]);
    diphoton_Photon_phi->Fill(PhotonPhi[0]);
    if (fabs(PhotonEtaSc[0]) < 1.442) {
      diphoton_Photon_hovere_eb->Fill(PhotonHoverE[0]);
      diphoton_Photon_sigmaietaieta_eb->Fill(PhotonSigmaIetaIeta[0]);
      diphoton_Photon_trksumptsolidconedr03_eb->Fill(PhotonTrkSumPtSolidConeDR03[0]);
      diphoton_Photon_e1x5e5x5_eb->Fill(PhotonE1x5E5x5[0]);
      diphoton_Photon_e2x5e5x5_eb->Fill(PhotonE2x5E5x5[0]);
    }
    if (fabs(PhotonEtaSc[0]) > 1.566 && fabs(PhotonEtaSc[0]) < 2.5) {
      diphoton_Photon_hovere_ee->Fill(PhotonHoverE[0]);
      diphoton_Photon_sigmaietaieta_ee->Fill(PhotonSigmaIetaIeta[0]);
      diphoton_Photon_trksumptsolidconedr03_ee->Fill(PhotonTrkSumPtSolidConeDR03[0]);
      diphoton_Photon_e1x5e5x5_ee->Fill(PhotonE1x5E5x5[0]);
      diphoton_Photon_e2x5e5x5_ee->Fill(PhotonE2x5E5x5[0]);
    }
    diphoton_Photon_energy->Fill(PhotonEnergy[1]);
    diphoton_Photon_pt->Fill(PhotonPt[1]);
    diphoton_Photon_et->Fill(PhotonEt[1]);
    diphoton_Photon_eta->Fill(PhotonEta[1]);
    diphoton_Photon_etasc->Fill(PhotonEtaSc[1]);
    diphoton_Photon_phi->Fill(PhotonPhi[1]);
    if (fabs(PhotonEtaSc[1]) < 1.4442) {
      diphoton_Photon_hovere_eb->Fill(PhotonHoverE[1]);
      diphoton_Photon_sigmaietaieta_eb->Fill(PhotonSigmaIetaIeta[1]);
      diphoton_Photon_trksumptsolidconedr03_eb->Fill(PhotonTrkSumPtSolidConeDR03[1]);
      diphoton_Photon_e1x5e5x5_eb->Fill(PhotonE1x5E5x5[1]);
      diphoton_Photon_e2x5e5x5_eb->Fill(PhotonE2x5E5x5[1]);
    }
    if (fabs(PhotonEtaSc[1]) > 1.566 && fabs(PhotonEtaSc[1]) < 2.5) {
      diphoton_Photon_hovere_ee->Fill(PhotonHoverE[1]);
      diphoton_Photon_sigmaietaieta_ee->Fill(PhotonSigmaIetaIeta[1]);
      diphoton_Photon_trksumptsolidconedr03_ee->Fill(PhotonTrkSumPtSolidConeDR03[1]);
      diphoton_Photon_e1x5e5x5_ee->Fill(PhotonE1x5E5x5[1]);
      diphoton_Photon_e2x5e5x5_ee->Fill(PhotonE2x5E5x5[1]);
    }
    diphoton_deltaPhiPhoton1Photon2->Fill(deltaPhi(PhotonPhi[0], PhotonPhi[1]));
    diphoton_deltaEtaPhoton1Photon2->Fill(PhotonEta[0] - PhotonEta[1]);
    diphoton_deltaRPhoton1Photon2->Fill(deltaR(PhotonEta[0], PhotonPhi[0], PhotonEta[1], PhotonPhi[1]));
    diphoton_invMassPhoton1Photon2->Fill(
        sqrt(2 * PhotonPt[0] * PhotonPt[1] * (cosh(PhotonEta[0] - PhotonEta[1]) - cos(PhotonPhi[0] - PhotonPhi[1]))));
    diphoton_PhotonMulti->Fill(diphoton_countPhoton_);
  }
}

void ExoticaDQM_miniAOD::analyzeMonoJets(const Event& iEvent) {
  for (unsigned int icoll = 0; icoll < MonoJetPFJetCollection_.size(); ++icoll) {
    monojet_countPFJet_ = 0;
    bool ValidMonoJetPFJets = iEvent.getByToken(MonoJetPFJetToken_[icoll], MonoJetpfJetCollection_);
    if (!ValidMonoJetPFJets)
      continue;
    MonoJetpfjets = *MonoJetpfJetCollection_;
    for (int i = 0; i < 2; i++) {
      PFJetPx[i] = 0.;
      PFJetPy[i] = 0.;
      PFJetPt[i] = 0.;
      PFJetEta[i] = 0.;
      PFJetPhi[i] = 0.;
      PFJetRapidity[i] = 0.;
      PFJetMass[i] = 0.;
      PFJetNHEF[i] = 0.;
      PFJetCHEF[i] = 0.;
      PFJetNEMF[i] = 0.;
      PFJetCEMF[i] = 0.;
    }
    std::vector<pat::Jet>::const_iterator MonoJetpfjet_ = MonoJetpfjets.begin();
    for (; MonoJetpfjet_ != MonoJetpfjets.end(); ++MonoJetpfjet_) {
      if (!MonoJetpfjet_->isPFJet()) { continue; }
      if (MonoJetpfjet_->pt() > PFJetPt[0]) {
        PFJetPt[0] = MonoJetpfjet_->pt();
        PFJetPx[0] = MonoJetpfjet_->px();
        PFJetPy[0] = MonoJetpfjet_->py();
        PFJetEta[0] = MonoJetpfjet_->eta();
        PFJetPhi[0] = MonoJetpfjet_->phi();
        PFJetRapidity[0] = MonoJetpfjet_->rapidity();
        PFJetMass[0] = MonoJetpfjet_->mass();
        PFJetNHEF[0] = MonoJetpfjet_->neutralHadronEnergyFraction();
        PFJetCHEF[0] = MonoJetpfjet_->chargedHadronEnergyFraction();
        PFJetNEMF[0] = MonoJetpfjet_->neutralEmEnergyFraction();
        PFJetCEMF[0] = MonoJetpfjet_->chargedEmEnergyFraction();
      } else {
      }
      if (MonoJetpfjet_->pt() > monojet_PFJet_pt_cut_)
        monojet_countPFJet_++;
    }
    if (PFJetPt[0] > monojet_PFJet_pt_cut_ && PFMETEt > monojet_PFJet_met_cut_) {
      monojet_PFJet_pt[icoll]->Fill(PFJetPt[0]);
      monojet_PFJet_eta[icoll]->Fill(PFJetEta[0]);
      monojet_PFJet_phi[icoll]->Fill(PFJetPhi[0]);
      monojet_PFMet[icoll]->Fill(PFMETEt);
      monojet_PFMet_phi[icoll]->Fill(PFMETPhi);
      monojet_PFJetPtOverPFMet[icoll]->Fill(PFJetPt[0] / PFMETEt);
      monojet_deltaPhiPFJetPFMet[icoll]->Fill(deltaPhi(PFJetPhi[0], PFMETPhi));
      monojet_PFchef[icoll]->Fill(PFJetCHEF[0]);
      monojet_PFnhef[icoll]->Fill(PFJetNHEF[0]);
      monojet_PFcemf[icoll]->Fill(PFJetCEMF[0]);
      monojet_PFnemf[icoll]->Fill(PFJetNEMF[0]);
      monojet_PFJetMulti[icoll]->Fill(monojet_countPFJet_);
    }
  }
}

void ExoticaDQM_miniAOD::analyzeMonoMuons(const Event& iEvent) {
  if (MuonPt[0] > monomuon_Muon_pt_cut_ && PFMETEt > monomuon_Muon_met_cut_) {
    monomuon_Muon_pt->Fill(MuonPt[0]);
    monomuon_Muon_eta->Fill(MuonEta[0]);
    monomuon_Muon_phi->Fill(MuonPhi[0]);
    monomuon_Charge->Fill(MuonCharge[0]);
    monomuon_PFMet->Fill(PFMETEt);
    monomuon_PFMet_phi->Fill(PFMETPhi);
    monomuon_MuonPtOverPFMet->Fill(MuonPt[0] / PFMETEt);
    monomuon_deltaPhiMuonPFMet->Fill(deltaPhi(MuonPhi[0], PFMETPhi));
    monomuon_TransverseMass->Fill(sqrt(2 * MuonPt[0] * PFMETEt * (1 - cos(deltaPhi(MuonPhi[0], PFMETPhi)))));
    monomuon_MuonMulti->Fill(monomuon_countMuon_);
  }
}

void ExoticaDQM_miniAOD::analyzeMonoElectrons(const Event& iEvent) {
  if (ElectronPt[0] > monoelectron_Electron_pt_cut_ && PFMETEt > monoelectron_Electron_met_cut_) {
    monoelectron_Electron_pt->Fill(ElectronPt[0]);
    monoelectron_Electron_eta->Fill(ElectronEta[0]);
    monoelectron_Electron_phi->Fill(ElectronPhi[0]);
    monoelectron_Charge->Fill(ElectronCharge[0]);
    monoelectron_PFMet->Fill(PFMETEt);
    monoelectron_PFMet_phi->Fill(PFMETPhi);
    monoelectron_ElectronPtOverPFMet->Fill(ElectronPt[0] / PFMETEt);
    monoelectron_deltaPhiElectronPFMet->Fill(deltaPhi(ElectronPhi[0], PFMETPhi));
    monoelectron_TransverseMass->Fill(
        sqrt(2 * ElectronPt[0] * PFMETEt * (1 - cos(deltaPhi(ElectronPhi[0], PFMETPhi)))));
    monoelectron_ElectronMulti->Fill(monoelectron_countElectron_);
  }
}

