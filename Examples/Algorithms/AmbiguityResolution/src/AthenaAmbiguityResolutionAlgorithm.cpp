// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolutionAlgorithm.hpp"
#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include <nlohmann/json.hpp>
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iterator>
#include <map>

namespace {

std::map<unsigned int,Acts::AthenaAmbiguityResolution::VolumeConfig>
readVolumeMap(std::string volumeFile) {
  std::ifstream file(volumeFile);
  if (!file.is_open()) {
      std::cerr << "Error opening file: " << volumeFile << std::endl;
      return {};
  }

  std::map<unsigned int,Acts::AthenaAmbiguityResolution::VolumeConfig> volumeMap;

  nlohmann::json j;

  file >> j;

  for (auto& [key, value] : j.items()) {

    unsigned int volumeId = std::stoi(key);

    int hitsScoreWeight = value["hitsScoreWeight"];
    int holesScoreWeight = value["holesScoreWeight"];
    int outliersScoreWeight = value["outliersScoreWeight"];
    int otherScoreWeight = value["otherScoreWeight"];

    std::size_t  minHits = value["minHits"];
    std::size_t  maxHits = value["maxHits"];
    std::size_t  maxHoles = value["maxHoles"];
    std::size_t  maxOutliers = value["maxOutliers"];
    std::size_t  maxUnused = value["maxUnused"];
    std::size_t  maxSharedHits = value["maxSharedHits"];

    bool sharedHitsFlag = value["sharedHitsFlag"];
    std::size_t detectorId = value["detectorId"];

    std::vector<double> goodHits = value["goodHits"];
    std::vector<double> fakeHits = value["fakeHits"];
    std::vector<double> goodHoles = value["goodHoles"];
    std::vector<double> fakeHoles = value["fakeHoles"];

    auto volumeConfig = Acts::AthenaAmbiguityResolution::VolumeConfig(hitsScoreWeight, holesScoreWeight, outliersScoreWeight,
                  otherScoreWeight, minHits, maxHits, maxHoles, maxOutliers, maxUnused, maxSharedHits, sharedHitsFlag, detectorId,
                  goodHits, fakeHits, goodHoles, fakeHoles);

    volumeMap[volumeId] = volumeConfig;
  }

  return volumeMap;
}


Acts::AthenaAmbiguityResolution::Config transformConfig(
    const ActsExamples::AthenaAmbiguityResolutionAlgorithm::Config& cfg, std::string volumeFile) {
  Acts::AthenaAmbiguityResolution::Config result;

  std::cout << "Volume File is "<<volumeFile << std::endl;

  result.volumeFile = volumeFile;
  result.volumeMap = readVolumeMap(result.volumeFile);
  result.minScore = cfg.minScore;
  result.minScoreSharedTracks = cfg.minScoreSharedTracks;
  result.maxSharedTracksPerMeasurement = cfg.maxSharedTracksPerMeasurement;
  result.maxShared = cfg.maxShared;
  result.pTMin = cfg.pTMin;
  result.pTMax = cfg.pTMax;
  result.phiMin = cfg.phiMin;
  result.phiMax = cfg.phiMax;
  result.etaMin = cfg.etaMin;
  result.etaMax = cfg.etaMax;
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
      m_core(transformConfig(cfg,m_cfg.volumeFile), logger().clone()) {
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

