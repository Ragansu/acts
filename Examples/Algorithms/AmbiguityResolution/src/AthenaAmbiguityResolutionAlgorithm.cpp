// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolutionAlgorithm.hpp"
#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iterator>
#include <map>

namespace {
Acts::AthenaAmbiguityResolution::Config transformConfig(
    const ActsExamples::AthenaAmbiguityResolutionAlgorithm::Config& cfg) {
  Acts::AthenaAmbiguityResolution::Config result;
    std::map<unsigned int,Acts::AthenaAmbiguityResolution::VolumeConfig> volumeMap = {

      {8,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 1
      {9,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 2 (barrel)
      {10,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 3

      {13,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 1 
      {14,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 2
      {15,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 3
      {16,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 4 (barrel)
      {18,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 5
      {19,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 6
      {20,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 7

      {22,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 1
      {23,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 2 (barrel)
      {24,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 3

      {25,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // HGTD 1
      {2,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // HGTD 2
    };
  std::cout<<"Volume map size: "<<cfg.volumeMap.size()<<std::endl;
  result.minScore = cfg.minScore;
  
  result.minScoreSharedTracks = cfg.minScoreSharedTracks;
  result.maxSharedTracksPerMeasurement = cfg.maxSharedTracksPerMeasurement;
  result.maxShared = cfg.maxShared;
  result.volumeMap = volumeMap;  // temporary hard-coded volume map
  return result;
}

std::size_t sourceLinkHash(const Acts::SourceLink& a) {
  return static_cast<std::size_t>(
      a.get<ActsExamples::IndexSourceLink>().index());
}

bool sourceLinkEquality(const Acts::SourceLink& a, const Acts::SourceLink& b) {
  return a.get<ActsExamples::IndexSourceLink>().index() ==
         b.get<ActsExamples::IndexSourceLink>().index();
}

}  // namespace

ActsExamples::AthenaAmbiguityResolutionAlgorithm::
    AthenaAmbiguityResolutionAlgorithm(
        ActsExamples::AthenaAmbiguityResolutionAlgorithm::Config cfg,
        Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("AthenaAmbiguityResolutionAlgorithm", lvl),
      m_cfg(std::move(cfg)),
      m_core(transformConfig(cfg), logger().clone()) {
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
  const auto& tracks = m_inputTracks(ctx);  // Read input data
  ACTS_VERBOSE("Number of input tracks: " << tracks.size());

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTracks;
  measurementsPerTracks = m_core.computeInitialState(tracks, &sourceLinkHash,
                              &sourceLinkEquality);

  std::vector<int> goodTracks =  m_core.solveAmbiguity(tracks, measurementsPerTracks);
  // Prepare the output track collection from the IDs
  TrackContainer solvedTracks{std::make_shared<Acts::VectorTrackContainer>(),
                            std::make_shared<Acts::VectorMultiTrajectory>()};
  solvedTracks.ensureDynamicColumns(tracks);
  for (auto iTrack : goodTracks) {
    auto destProxy = solvedTracks.getTrack(solvedTracks.addTrack());
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFrom(srcProxy, false);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  ActsExamples::ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(solvedTracks.container())),
      tracks.trackStateContainerHolder()};

  m_outputTracks(ctx, std::move(outputTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}

