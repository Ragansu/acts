// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolutionAlgorithm.hpp"

#include "ActsExamples/Framework/ProcessCode.hpp"

#include <iterator>
#include <map>

ActsExamples::AthenaAmbiguityResolutionAlgorithm::AthenaAmbiguityResolutionAlgorithm(
    ActsExamples::AthenaAmbiguityResolutionAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::AthenaAmbiguityResolution("AthenaAmbiguityResolutionAlgorithm",
                                           lvl),
      m_cfg(std::move(cfg)){
  
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode ActsExamples::AthenaAmbiguityResolutionAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& tracks = m_inputTracks(ctx);
  // Associate measurement to their respective tracks
  std::map<std::size_t, Counter> counterMap;
  std::vector<int> score = simpleScore(tracks, counterMap);
  // Select the ID of the track we want to keep
  std::vector<std::size_t> goodTracks = solveAmbiguity(tracks,score,counterMap);
  // Prepare the output track collection from the IDs
  auto outputTracks = prepareOutputTrack(tracks, goodTracks);
  m_outputTracks(ctx, std::move(outputTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}