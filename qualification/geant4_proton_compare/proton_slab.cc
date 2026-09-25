#include "G4Box.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Proton.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "G4UserLimits.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4hIonisation.hh"
#include "Randomize.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
double slab_cm = 0.01;
constexpr double source_half_width_MeV = 0.005;
constexpr std::array<double, 3> energies_MeV{10.0, 100.0, 499.99};

class Detector final : public G4VUserDetectorConstruction {
 public:
  explicit Detector(double max_step_mm) : max_step_mm_(max_step_mm) {}
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    auto* al = nist->FindOrBuildMaterial("G4_Al");
    auto* cu = nist->FindOrBuildMaterial("G4_Cu");
    auto* world = new G4LogicalVolume(new G4Box("world", 10*cm, 10*cm, 10*cm), vacuum, "world");
    auto* world_pv = new G4PVPlacement(nullptr, {}, world, "world", nullptr, false, 0);
    auto* al_lv = new G4LogicalVolume(new G4Box("al", slab_cm*cm/2, 10*cm, 10*cm), al, "al");
    auto* cu_lv = new G4LogicalVolume(new G4Box("cu", slab_cm*cm/2, 10*cm, 10*cm), cu, "cu");
    new G4PVPlacement(nullptr, G4ThreeVector(slab_cm*cm/2, 0, 0), al_lv, "al", world, false, 0);
    new G4PVPlacement(nullptr, G4ThreeVector(3*slab_cm*cm/2, 0, 0), cu_lv, "cu", world, false, 0);
    al_lv->SetUserLimits(new G4UserLimits(max_step_mm_*mm));
    cu_lv->SetUserLimits(new G4UserLimits(max_step_mm_*mm));
    return world_pv;
  }
 private:
  double max_step_mm_;
};

class IonisationOnly final : public G4VUserPhysicsList {
 public:
  IonisationOnly() {
    defaultCutValue = 1*m;
    G4EmParameters::Instance()->SetLossFluctuations(false);
    G4EmParameters::Instance()->SetStepFunction(0.01, 0.001*mm);
  }
  void ConstructParticle() override {
    G4Proton::ProtonDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4Gamma::GammaDefinition();
  }
  void ConstructProcess() override {
    AddTransportation();
    auto* manager = G4Proton::ProtonDefinition()->GetProcessManager();
    ion_ = new G4hIonisation();
    manager->AddProcess(ion_, -1, 2, 2);
    manager->AddDiscreteProcess(new G4StepLimiter());
  }
  void SetCuts() override { SetCutsWithDefault(); }
  G4hIonisation* Ionisation() const { return ion_; }
 private:
  G4hIonisation* ion_ = nullptr;
};

class Primary final : public G4VUserPrimaryGeneratorAction {
 public:
  Primary() : gun_(new G4ParticleGun(1)) { gun_->SetParticleDefinition(G4Proton::ProtonDefinition()); }
  ~Primary() override { delete gun_; }
  void SetEnergy(double e) { energy_MeV_ = e; }
  void GeneratePrimaries(G4Event* event) override {
    const double mu = 1/std::sqrt(3.0);
    const double sign = G4UniformRand() < 0.5 ? -1.0 : 1.0;
    gun_->SetParticlePosition(G4ThreeVector(G4UniformRand()*slab_cm*cm, 0, 0));
    gun_->SetParticleMomentumDirection(G4ThreeVector(sign*mu, std::sqrt(2.0/3.0), 0));
    gun_->SetParticleEnergy((energy_MeV_ + (2*G4UniformRand()-1)*source_half_width_MeV)*MeV);
    gun_->GeneratePrimaryVertex(event);
  }
 private:
  G4ParticleGun* gun_;
  double energy_MeV_ = 10;
};

class Scores final {
 public:
  void Reset() { sum_ = {0,0}; sumsq_ = {0,0}; count_ = 0; }
  void Begin() { event_ = {0,0}; }
  void Step(const G4Step* step) {
    const auto* volume = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!volume) return;
    const auto& name = volume->GetName();
    if (name == "al") event_[0] += step->GetTotalEnergyDeposit()/MeV;
    if (name == "cu") event_[1] += step->GetTotalEnergyDeposit()/MeV;
  }
  void End() {
    ++count_;
    for (int i=0; i<2; ++i) { sum_[i] += event_[i]; sumsq_[i] += event_[i]*event_[i]; }
  }
  std::array<double,2> Mean() const { return {sum_[0]/count_, sum_[1]/count_}; }
  std::array<double,2> StandardError() const {
    std::array<double,2> out{};
    for (int i=0; i<2; ++i) {
      const double mean = sum_[i]/count_;
      out[i] = std::sqrt(std::max(0.0, (sumsq_[i]/count_ - mean*mean)/(count_-1)));
    }
    return out;
  }
 private:
  std::array<double,2> event_{}, sum_{}, sumsq_{};
  int count_ = 0;
};

class EventScores final : public G4UserEventAction {
 public:
  explicit EventScores(Scores& scores) : scores_(scores) {}
  void BeginOfEventAction(const G4Event*) override { scores_.Begin(); }
  void EndOfEventAction(const G4Event*) override { scores_.End(); }
 private:
  Scores& scores_;
};

class StepScores final : public G4UserSteppingAction {
 public:
  explicit StepScores(Scores& scores) : scores_(scores) {}
  void UserSteppingAction(const G4Step* step) override { scores_.Step(step); }
 private:
  Scores& scores_;
};
}

int main(int argc, char** argv) {
  if (argc != 4 && argc != 5 && argc != 7) {
    std::cerr << "usage: radiant_proton_compare <stopping.csv> <scores.csv> <histories-per-energy> [max-step-mm [single-energy-MeV slab-thickness-cm]]\n";
    return 2;
  }
  const int histories = std::stoi(argv[3]);
  const double max_step_mm = argc >= 5 ? std::stod(argv[4]) : 0.005;
  std::vector<double> run_energies(energies_MeV.begin(), energies_MeV.end());
  if (argc == 7) {
    run_energies = {std::stod(argv[5])};
    slab_cm = std::stod(argv[6]);
  }
  if (histories < 2 || !(max_step_mm > 0) || !(slab_cm > 0) ||
      std::any_of(run_energies.begin(), run_energies.end(), [](double e){ return e - source_half_width_MeV <= 1 || e + source_half_width_MeV > 500; })) return 2;
  CLHEP::HepRandom::setTheSeed(731293);
  auto* manager = new G4RunManager();
  manager->SetUserInitialization(new Detector(max_step_mm));
  auto* physics = new IonisationOnly();
  manager->SetUserInitialization(physics);
  auto* primary = new Primary();
  Scores scores;
  manager->SetUserAction(primary);
  manager->SetUserAction(new EventScores(scores));
  manager->SetUserAction(new StepScores(scores));
  manager->Initialize();
  std::cerr << "Geant4 initialized; warming one event\n";
  primary->SetEnergy(100.0);
  manager->BeamOn(1);
  std::cerr << "Geant4 event loop initialized; computing dE/dx\n";

  auto* nist = G4NistManager::Instance();
  auto* al = nist->FindOrBuildMaterial("G4_Al");
  auto* cu = nist->FindOrBuildMaterial("G4_Cu");
  auto* cuts = G4ProductionCutsTable::GetProductionCutsTable();
  const G4MaterialCutsCouple* al_couple = nullptr;
  const G4MaterialCutsCouple* cu_couple = nullptr;
  for (std::size_t i=0; i<cuts->GetTableSize(); ++i) {
    const auto* couple = cuts->GetMaterialCutsCouple(static_cast<int>(i));
    if (couple->GetMaterial() == al) al_couple = couple;
    if (couple->GetMaterial() == cu) cu_couple = couple;
  }
  if (!al_couple || !cu_couple || !physics->Ionisation()) {
    throw std::runtime_error("Missing Geant4 material/cut couple or ionisation process");
  }
  std::ofstream table(argv[1]);
  table << "energy_MeV,Al_MeV_cm,Cu_MeV_cm\n" << std::setprecision(15);
  for (int i=0; i<=300; ++i) {
    const double energy = std::exp(std::log(1.0) + i*std::log(500.0)/300);
    const double al_dedx = physics->Ionisation()->GetDEDX(energy*MeV, al_couple)/(MeV/cm);
    const double cu_dedx = physics->Ionisation()->GetDEDX(energy*MeV, cu_couple)/(MeV/cm);
    table << energy << ',' << al_dedx << ',' << cu_dedx << '\n';
  }
  table.close();
  std::ofstream out(argv[2]);
  out << "energy_MeV,histories,Al_density_g_cm3,Cu_density_g_cm3,Al_edep_MeV,Al_se_MeV,Cu_edep_MeV,Cu_se_MeV,max_step_mm,slab_cm\n";
  out << std::setprecision(15);
  for (const double energy : run_energies) {
    scores.Reset();
    primary->SetEnergy(energy);
    manager->BeamOn(histories);
    const auto mean = scores.Mean();
    const auto se = scores.StandardError();
    out << energy << ',' << histories << ',' << al->GetDensity()/(g/cm3) << ',' << cu->GetDensity()/(g/cm3)
        << ',' << mean[0] << ',' << se[0] << ',' << mean[1] << ',' << se[1] << ',' << max_step_mm << ',' << slab_cm << '\n';
    std::cout << "E=" << energy << " MeV Al=" << mean[0] << " Cu=" << mean[1] << "\n";
  }
  delete manager;
  return 0;
}
