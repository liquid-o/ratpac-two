#ifndef __RAT_GeoFiberSensitiveDetectorHit__
#define __RAT_GeoFiberSensitiveDetectorHit__

#include <G4Allocator.hh>
#include <G4LogicalVolume.hh>
#include <G4RotationMatrix.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VHit.hh>

namespace RAT {

class GeoFiberSensitiveDetectorHit : public G4VHit {
 public:
  GeoFiberSensitiveDetectorHit(G4int i, G4double t, G4ThreeVector p);
  virtual ~GeoFiberSensitiveDetectorHit();
  GeoFiberSensitiveDetectorHit(const GeoFiberSensitiveDetectorHit &right);
  const GeoFiberSensitiveDetectorHit &operator=(const GeoFiberSensitiveDetectorHit &right);
  int operator==(const GeoFiberSensitiveDetectorHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Draw();
  void Print();

 private:
  G4int id;
  G4double time;
  G4ThreeVector pos;
  G4ThreeVector hit_pos;
  G4RotationMatrix rot;
  const G4LogicalVolume *pLogV;
  std::string proc;
  std::vector<G4int> readout_ch_nums = {0, 0};

 public:
  inline G4int GetID() const { return id; }
  inline std::vector<G4int> GetReadoutChNums() const { return readout_ch_nums; }
  inline void SetReadoutChNums(G4int val1, G4int val2) { readout_ch_nums[0] = val1; readout_ch_nums[1] = val2; }
  inline G4double GetTime() const { return time; }
  inline void SetTime(G4double val) { time = val; }
  inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
  inline G4ThreeVector GetPos() const { return pos; }
  inline void SetHitPos(G4ThreeVector xyz) { hit_pos = xyz; }
  inline G4ThreeVector GetHitPos() const { return hit_pos; }
  inline void SetRot(G4RotationMatrix rmat) { rot = rmat; }
  inline G4RotationMatrix GetRot() const { return rot; }
  inline void SetLogV(G4LogicalVolume *val) { pLogV = val; }
  inline const G4LogicalVolume *GetLogV() const { return pLogV; }
  inline void SetProc(std::string val) { proc = val; }
  inline const std::string GetProc() const { return proc; }
};

typedef G4THitsCollection<GeoFiberSensitiveDetectorHit> GeoFiberSensitiveDetectorHitsCollection;

extern G4Allocator<GeoFiberSensitiveDetectorHit> GeoFiberSensitiveDetectorHitAllocator;

inline void *GeoFiberSensitiveDetectorHit::operator new(size_t) {
  void *aHit;
  aHit = (void *)GeoFiberSensitiveDetectorHitAllocator.MallocSingle();
  return aHit;
}

inline void GeoFiberSensitiveDetectorHit::operator delete(void *aHit) {
  GeoFiberSensitiveDetectorHitAllocator.FreeSingle((GeoFiberSensitiveDetectorHit *)aHit);
}

}  // namespace RAT

#endif
