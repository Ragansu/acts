// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/AmbiguityResolution/GreedyAmbiguityResolutionAlgorithm.hpp"
#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolutionAlgorithm.hpp"


#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {

void addAmbiguityResolution(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::GreedyAmbiguityResolutionAlgorithm, mex,
      "GreedyAmbiguityResolutionAlgorithm", inputTracks, outputTracks,
      maximumSharedHits, maximumIterations, nMeasurementsMin);
}
void addAthenaAmbiguityResolution(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::AthenaAmbiguityResolutionAlgorithm, mex,
      "AthenaAmbiguityResolutionAlgorithm", inputTracks, outputTracks,
      minScore, minScoreSharedTracks, maxShared, maxSharedTracksPerMeasurement,
      phiMin, phiMax, etaMin, etaMax);
}

}  // namespace Acts::Python
