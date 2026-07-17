////////////////////////////////////////////////////////////////////
/// \class RAT::NoiseReductionNNProc
///
/// \brief Nearest-neighbor noise reduction for digitized SiPM channels
///
/// \details
/// Removes channels whose surviving fiber neighborhood is too sparse
/// to be consistent with a light signal. Operates at the *fiber* level:
/// at BeginOfRun, channels are paired into fibers (two channels per
/// fiber, one at each z-end, matched by (x,y) position) and each fiber
/// gets a neighbor list of other fibers within [d_min, d_min * ratio]
/// in the (x,y) plane.
///
/// Per-event cut (single pass, counts computed against the original
/// surviving set): a channel is kept iff the sum of surviving hits on
/// its fiber's neighbor fibers is >= nn_required. A neighbor fiber
/// with hits at both ends contributes 2 to the count; a single-hit
/// neighbor contributes 1. This is a pure nearest-neighbor cut in the
/// (x,y) plane.
///
/// Optional Rule A "two-sided keep" (legacy): when disable_rule_a = 0,
/// any channel on a fiber that has hits at BOTH z-ends is auto-kept
/// regardless of its neighbor count. This was the original behavior
/// but was found to inflate noise survival from accidental dark-count
/// two-ended coincidences. Default disable_rule_a = 1 (auto-keep off).
///
/// Parameters are loaded from NOISE_REDUCTION_NN.ratdb:
///   enabled, nn_required (default 3), max_nn_distance_ratio,
///   diagnostics, disable_rule_a (default 1).
///
/// Intended placement: directly after WaveformPrep in the processor
/// chain, which has already applied the voltage-threshold zero-suppress.
////////////////////////////////////////////////////////////////////
#ifndef __RAT_NoiseReductionNNProc__
#define __RAT_NoiseReductionNNProc__

#include <RAT/DB.hh>
#include <RAT/Processor.hh>

#include <string>
#include <vector>

namespace RAT {

namespace DS {
class Root;
class EV;
class Run;
}  // namespace DS

class NoiseReductionNNProc : public Processor {
 public:
  NoiseReductionNNProc();
  virtual ~NoiseReductionNNProc() {}

  virtual void BeginOfRun(DS::Run *run);
  virtual Processor::Result Event(DS::Root *ds, DS::EV *ev);

  virtual void SetI(std::string param, int value);
  virtual void SetD(std::string param, double value);

 protected:
  void LoadParameters();
  void BuildFibersAndNeighbors(DS::Run *run);

  DBLinkPtr fLink;

  // ratdb-configured
  int fEnabled;
  int fNNRequired;
  double fMaxNNDistanceRatio;
  int fDiagnostics;
  // When non-zero, skip the Rule A short-circuit so a fiber is kept
  // iff its Rule B neighbor count >= nn_required (pure NN cut in xy).
  // Both-end fibers still contribute 2 to OTHER fibers' counts.
  int fDisableRuleA;

  // Geometry caches, filled at BeginOfRun.
  // fFiberOfChannel[pmtID] -> fiber ID, or -1 if the channel was not
  // seen during the geometry scan.
  std::vector<int> fFiberOfChannel;
  // fFiberChannels[fiber ID] -> list of channel pmtIDs on that fiber
  // (size 1 for a single-ended fiber, 2 for a two-ended fiber).
  std::vector<std::vector<int>> fFiberChannels;
  // fFiberNeighbors[fiber ID] -> list of neighbor fiber IDs.
  std::vector<std::vector<int>> fFiberNeighbors;
  bool fNeighborMapBuilt;
};

}  // namespace RAT

#endif
