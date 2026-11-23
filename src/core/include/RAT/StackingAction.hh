// fsutanto@umich.edu
// Apr 15,2018
// the class is added to kill secondary tracks
// that are produced by neutron capture on 158Gd

#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

#include <vector>

struct OpPhotonRecord {
  int trackID, parentID;
  short creatorProc;   // 0: scint, 1: Cerenkov, 2: reemission, 3: other
  float x, y, z, t;    // mm, ns
  float dx, dy, dz;
  float px, py, pz;
  float energy;        // MeV
  float weight;
};

class StackingAction : public G4UserStackingAction {
 public:
  StackingAction();
  virtual ~StackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *);
  void PrepareNewEvent() override;
  const std::vector<OpPhotonRecord> &Photons() const { return fPhotons; }
  bool DumpEnabled() const { return fDumpOpticalPhotons; }

 private:
  bool fDumpOpticalPhotons = false;
  std::vector<OpPhotonRecord> fPhotons;
};

#endif
