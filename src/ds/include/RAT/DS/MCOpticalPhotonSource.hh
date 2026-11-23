#ifndef __RAT_DS_MCOpticalPhotonSource__
#define __RAT_DS_MCOpticalPhotonSource__

namespace RAT {
namespace DS {

/**
 * @struct DS::MCOpticalPhotonSource
 * Optical photon recorded (and killed) at creation, for hand-off to a fast
 * optical simulation. Transient: never written to the standard output.
 */
struct MCOpticalPhotonSource {
  enum Process { kScintillation = 0, kCerenkov = 1, kReemission = 2, kOther = 3 };

  int trackID = 0;
  int parentID = 0;
  short creatorProc = kOther;
  float x = 0, y = 0, z = 0;     // mm
  float t = 0;                   // ns, global time (includes scintillation delay)
  float dx = 0, dy = 0, dz = 0;  // unit momentum direction
  float px = 0, py = 0, pz = 0;  // polarization
  float energy = 0;              // MeV
  float weight = 1;
};

}  // namespace DS
}  // namespace RAT

#endif
