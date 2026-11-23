// fsutanto@umich.edu
// Apr 15,2018
// the class is added to kill secondary tracks
// that are produced by neutron capture on 158Gd

// Add option to kill photons upon production, recording their information

#include "RAT/StackingAction.hh"

#include <RAT/DB.hh>
#include <RAT/EventInfo.hh>
#include <RAT/Log.hh>
#include <RAT/TrackInfo.hh>

#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Gamma.hh"
#include "G4NeutrinoE.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

// constructor
StackingAction::StackingAction() {}

// destructor
StackingAction::~StackingAction() {}

void StackingAction::PrepareNewEvent() {
  fPhotons.clear();
  try {
    fDumpOpticalPhotons = RAT::DB::Get()->GetLink("MC")->GetI("dump_optical_photons") != 0;
  } catch (RAT::DBNotFoundError &) {
    fDumpOpticalPhotons = false;
  }
}

// Classify new tracks
G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track *track) {
  G4ParticleDefinition *par = track->GetDefinition();
  G4String nameParticle = par->GetParticleName();

#define debug_dicebox
#undef debug_dicebox
#ifdef debug_dicebox

  if (nameParticle != "opticalphoton" && track->GetParentID() == 1) {
    RAT::debug << "*******************************" << newline;
    RAT::debug << "      Classify new track       " << newline;
    RAT::debug << "*******************************" << newline;

    RAT::debug << "(name, parentID,trackID,stepID,status) : (" << par->GetParticleName() << "," << track->GetParentID()
               << "," << track->GetTrackID() << "," << track->GetCurrentStepNumber() << "," << track->GetTrackStatus()
               << ")" << newline;
  }

#endif

  // Optical photon dump: record, count, kill
  if (fDumpOpticalPhotons && par == G4OpticalPhoton::OpticalPhotonDefinition()) {
    // same creator-process logic as Gsim::PreUserTrackingAction
    std::string proc;
    if (const G4VProcess *cp = track->GetCreatorProcess()) proc = cp->GetProcessName();
    auto *ti = dynamic_cast<RAT::TrackInfo *>(track->GetUserInformation());
    if (ti && ti->GetCreatorProcess() != "") proc = ti->GetCreatorProcess();

    G4Event *event = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
    auto *ei = event ? dynamic_cast<RAT::EventInfo *>(event->GetUserInformation()) : nullptr;

    OpPhotonRecord r;
    if (proc.find("Scintillation") != std::string::npos) {
      r.creatorProc = 0;
      if (ei) ei->numScintPhoton++;
    } else if (proc.find("Cerenkov") != std::string::npos) {
      r.creatorProc = 1;
      if (ei) ei->numCerenkovPhoton++;
    } else if (proc.find("Reemission") != std::string::npos) {
      r.creatorProc = 2;
      if (ei) ei->numReemitPhoton++;
    } else {
      r.creatorProc = 3;
    }

    r.trackID = track->GetTrackID();
    r.parentID = track->GetParentID();
    const G4ThreeVector &pos = track->GetPosition();
    r.x = pos.x(); r.y = pos.y(); r.z = pos.z();
    r.t = track->GetGlobalTime();
    const G4ThreeVector &d = track->GetMomentumDirection();
    r.dx = d.x(); r.dy = d.y(); r.dz = d.z();
    const G4ThreeVector &pol = track->GetPolarization();
    r.px = pol.x(); r.py = pol.y(); r.pz = pol.z();
    r.energy = track->GetKineticEnergy();
    r.weight = track->GetWeight();
    fPhotons.push_back(r);
    return fKill;
  }

  // In SteppingAction.cc, if a step is neutron capture on 157Gd,
  // we postpone all of its secondaries except the 158Gd.
  // These secondaries are then send here (StackingAction.cc)
  // to be removed.
  // GetTrackStatus()==5 is the fPostponeToNextEvent

  if (track->GetTrackStatus() == 5 && (nameParticle == "gamma" || nameParticle == "e-")) {
    return fKill;
  } else {
    return fUrgent;
  }
}
