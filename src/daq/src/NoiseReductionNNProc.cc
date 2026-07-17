#include <TVector3.h>

#include <RAT/DS/EV.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/Log.hh>
#include <RAT/NoiseReductionNNProc.hh>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace RAT {

NoiseReductionNNProc::NoiseReductionNNProc()
    : Processor("noisereductionnn"),
      fEnabled(0),
      fNNRequired(0),
      fMaxNNDistanceRatio(1.45),
      fDiagnostics(0),
      fDisableRuleA(1),
      fNeighborMapBuilt(false) {
  LoadParameters();
}

void NoiseReductionNNProc::LoadParameters() {
  try {
    fLink = DB::Get()->GetLink("NOISE_REDUCTION_NN");
    fEnabled = fLink->GetI("enabled");
    fNNRequired = fLink->GetI("nn_required");
    fMaxNNDistanceRatio = fLink->GetD("max_nn_distance_ratio");
    fDiagnostics = fLink->GetI("diagnostics");
  } catch (DBNotFoundError &e) {
    Log::Die("NoiseReductionNNProc: Unable to find NOISE_REDUCTION_NN in ratdb.");
  }
  // Optional key — default 0 if absent so existing ratdb files keep working.
  try { fDisableRuleA = fLink->GetI("disable_rule_a"); } catch (DBNotFoundError &e) {}
}

void NoiseReductionNNProc::SetI(std::string param, int value) {
  if (param == "enabled") {
    fEnabled = value;
  } else if (param == "nn_required") {
    fNNRequired = value;
  } else if (param == "diagnostics") {
    fDiagnostics = value;
  } else if (param == "disable_rule_a") {
    fDisableRuleA = value;
  } else {
    throw Processor::ParamUnknown(param);
  }
}

void NoiseReductionNNProc::SetD(std::string param, double value) {
  if (param == "max_nn_distance_ratio") {
    if (value != fMaxNNDistanceRatio) fNeighborMapBuilt = false;
    fMaxNNDistanceRatio = value;
  } else {
    throw Processor::ParamUnknown(param);
  }
}

void NoiseReductionNNProc::BeginOfRun(DS::Run *run) {
  if (!fEnabled) return;
  BuildFibersAndNeighbors(run);
}

void NoiseReductionNNProc::BuildFibersAndNeighbors(DS::Run *run) {
  fFiberOfChannel.clear();
  fFiberChannels.clear();
  fFiberNeighbors.clear();
  fNeighborMapBuilt = false;

  DS::PMTInfo *pmtinfo = run->GetPMTInfo();
  if (!pmtinfo) {
    warn << "NoiseReductionNNProc: PMTInfo unavailable at BeginOfRun; skipping." << newline;
    return;
  }

  const int nPMT = pmtinfo->GetPMTCount();
  if (nPMT <= 0) {
    warn << "NoiseReductionNNProc: no PMTs found in PMTInfo." << newline;
    return;
  }

  // Cache positions and determine z-end split relative to mean(z).
  std::vector<TVector3> pos(nPMT);
  double mean_z = 0.0;
  for (int i = 0; i < nPMT; ++i) {
    pos[i] = pmtinfo->GetPosition(i);
    mean_z += pos[i].Z();
  }
  mean_z /= nPMT;

  std::vector<int> posIds, negIds;
  posIds.reserve(nPMT / 2);
  negIds.reserve(nPMT / 2);
  for (int i = 0; i < nPMT; ++i) {
    if (pos[i].Z() >= mean_z)
      posIds.push_back(i);
    else
      negIds.push_back(i);
  }

  // Pair each +z channel with its opposite-end partner at the same (x,y).
  // 10 μm tolerance easily covers float<->double rounding on positions that
  // are constructed identically per fiber.
  const double pair_tol = 1e-2;
  const double pair_tol2 = pair_tol * pair_tol;
  std::vector<int> partner(nPMT, -1);
  for (int pi : posIds) {
    double bestD2 = std::numeric_limits<double>::infinity();
    int bestJ = -1;
    for (int nj : negIds) {
      double dx = pos[pi].X() - pos[nj].X();
      double dy = pos[pi].Y() - pos[nj].Y();
      double d2 = dx * dx + dy * dy;
      if (d2 < bestD2) {
        bestD2 = d2;
        bestJ = nj;
      }
    }
    if (bestJ >= 0 && bestD2 <= pair_tol2 && partner[bestJ] == -1) {
      partner[pi] = bestJ;
      partner[bestJ] = pi;
    }
  }

  // Assign fiber IDs. A single-ended channel gets its own one-channel fiber.
  fFiberOfChannel.assign(nPMT, -1);
  for (int i = 0; i < nPMT; ++i) {
    if (fFiberOfChannel[i] >= 0) continue;
    int fid = static_cast<int>(fFiberChannels.size());
    fFiberChannels.push_back({i});
    fFiberOfChannel[i] = fid;
    if (partner[i] >= 0) {
      fFiberOfChannel[partner[i]] = fid;
      fFiberChannels.back().push_back(partner[i]);
    }
  }

  const int nFibers = static_cast<int>(fFiberChannels.size());

  // Fiber (x,y): first channel of each fiber.
  std::vector<double> fx(nFibers), fy(nFibers);
  for (int fid = 0; fid < nFibers; ++fid) {
    int c0 = fFiberChannels[fid][0];
    fx[fid] = pos[c0].X();
    fy[fid] = pos[c0].Y();
  }

  // Per-fiber nearest-neighbor distance in (x,y).
  std::vector<double> dmin(nFibers, std::numeric_limits<double>::infinity());
  for (int i = 0; i < nFibers; ++i) {
    for (int j = 0; j < nFibers; ++j) {
      if (j == i) continue;
      double dx = fx[i] - fx[j];
      double dy = fy[i] - fy[j];
      double d = std::sqrt(dx * dx + dy * dy);
      if (d > 0.0 && d < dmin[i]) dmin[i] = d;
    }
  }

  // Build fiber neighbor list: fibers with d in [d_min, d_min * ratio].
  fFiberNeighbors.assign(nFibers, std::vector<int>());
  const double tol = 1e-6;
  for (int i = 0; i < nFibers; ++i) {
    if (!std::isfinite(dmin[i])) continue;
    const double d_upper = dmin[i] * fMaxNNDistanceRatio;
    for (int j = 0; j < nFibers; ++j) {
      if (j == i) continue;
      double dx = fx[i] - fx[j];
      double dy = fy[i] - fy[j];
      double d = std::sqrt(dx * dx + dy * dy);
      if (d >= dmin[i] - tol && d <= d_upper + tol) fFiberNeighbors[i].push_back(j);
    }
  }

  fNeighborMapBuilt = true;

  if (fDiagnostics) {
    int n_single = 0, n_pair = 0;
    for (auto const &chs : fFiberChannels) {
      if (chs.size() == 1)
        n_single++;
      else if (chs.size() == 2)
        n_pair++;
    }
    long total_nb = 0;
    int n_with_nb = 0;
    for (auto const &v : fFiberNeighbors) {
      if (!v.empty()) {
        total_nb += v.size();
        n_with_nb++;
      }
    }
    const double mean_nb = (n_with_nb > 0) ? static_cast<double>(total_nb) / n_with_nb : 0.0;
    info << "NoiseReductionNNProc: " << nPMT << " PMTs -> " << nFibers << " fibers (" << n_pair << " two-ended, "
         << n_single << " single-ended); " << n_with_nb << " fibers with >=1 neighbor; mean neighbor-fiber count = "
         << mean_nb << "; nn_required=" << fNNRequired << "; max_nn_distance_ratio=" << fMaxNNDistanceRatio << newline;
  }
}

Processor::Result NoiseReductionNNProc::Event(DS::Root * /*ds*/, DS::EV *ev) {
  if (!fEnabled) return Processor::Result::OK;
  if (!fNeighborMapBuilt) return Processor::Result::OK;
  if (fNNRequired <= 0) return Processor::Result::OK;

  const std::vector<int> ids = ev->GetAllDigitPMTIDs();
  if (ids.empty()) return Processor::Result::OK;

  // Aggregate surviving channel hits per fiber.
  std::map<int, int> hits_per_fiber;
  for (int id : ids) {
    if (id < 0 || id >= static_cast<int>(fFiberOfChannel.size())) continue;
    int fid = fFiberOfChannel[id];
    if (fid >= 0) hits_per_fiber[fid]++;
  }

  // Decide which channels to erase, using the single-pass rules.
  std::vector<int> toErase;
  toErase.reserve(ids.size());
  for (int id : ids) {
    if (id < 0 || id >= static_cast<int>(fFiberOfChannel.size())) {
      toErase.push_back(id);
      continue;
    }
    int fid = fFiberOfChannel[id];
    if (fid < 0) {
      toErase.push_back(id);
      continue;
    }
    // Rule A: a fiber hit at both ends is always kept (unless gated off).
    // When fDisableRuleA != 0, both-end fibers still contribute 2 to neighbor
    // counts (Rule B), but are no longer auto-kept here — they must satisfy
    // Rule B on their own to survive.
    if (!fDisableRuleA && hits_per_fiber[fid] >= 2) continue;
    // Rule B: single-hit fiber — count surviving hits on neighbor fibers.
    int count = 0;
    for (int nfid : fFiberNeighbors[fid]) {
      auto it = hits_per_fiber.find(nfid);
      if (it != hits_per_fiber.end()) count += it->second;
    }
    if (count < fNNRequired) toErase.push_back(id);
  }

  for (int id : toErase) ev->EraseDigitPMT(id);

  if (fDiagnostics) {
    debug << "NoiseReductionNNProc: before=" << ids.size() << " erased=" << toErase.size()
          << " after=" << (ids.size() - toErase.size()) << newline;
  }

  return Processor::Result::OK;
}

}  // namespace RAT

