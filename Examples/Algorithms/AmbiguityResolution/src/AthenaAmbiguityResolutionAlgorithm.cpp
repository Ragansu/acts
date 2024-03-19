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

    int minHits = value["minHits"];
    int maxHoles = value["maxHoles"];
    int maxOutliers = value["maxOutliers"];
    int maxUnused = value["maxUnused"];
    int maxSharedHits = value["maxSharedHits"];

    bool sharedHitsFlag = value["sharedHitsFlag"];
    std::size_t detectorId = value["detectorId"];


    volumeMap[volumeId] = {
      hitsScoreWeight, 
      holesScoreWeight, 
      outliersScoreWeight, 
      otherScoreWeight,
      minHits, 
      maxHoles, 
      maxOutliers, 
      maxUnused, 
      maxSharedHits, 
      sharedHitsFlag, 
      detectorId
    };
  }

  return volumeMap;
}


Acts::AthenaAmbiguityResolution::Config transformConfig(
    const ActsExamples::AthenaAmbiguityResolutionAlgorithm::Config& cfg) {
  Acts::AthenaAmbiguityResolution::Config result;
  result.volumeMap = readVolumeMap(cfg.volumeFile);
  result.minScore = cfg.minScore;
  result.minScoreSharedTracks = cfg.minScoreSharedTracks;
  result.maxSharedTracksPerMeasurement = cfg.maxSharedTracksPerMeasurement;
  result.maxShared = cfg.maxShared;
  result.phiMin = cfg.phiMin;
  result.phiMax = cfg.phiMax;
  result.etaMin = cfg.etaMin;
  result.etaMax = cfg.etaMax;
  return result;
}


// std::map<unsigned int,Acts::AthenaAmbiguityResolution::VolumeConfig>
// readVolumeMap(std::string volumeFile) {
//   std::ifstream file(volumeFile);
//   if (!file.is_open()) {
//       std::cerr << "Error opening file: " << volumeFile << std::endl;
//       return {};
//   }

//   std::map<unsigned int,Acts::AthenaAmbiguityResolution::VolumeConfig> volumeMap;

//   std::string line;
//   while (std::getline(file, line)) {
//     std::istringstream iss(line);
//     std::string key_str, value_str;
//     if (std::getline(iss, key_str, ':') && std::getline(iss, value_str)) {
//       // Remove leading/trailing whitespace
//       key_str.erase(0, key_str.find_first_not_of(" \t"));
//       key_str.erase(key_str.find_last_not_of(" \t") + 1);
//       value_str.erase(0, value_str.find_first_not_of(" \t"));
//       value_str.erase(value_str.find_last_not_of(" \t") + 1);

//       unsigned int volumeId = std::stoi(key_str);

//       int hitsScoreWeight;
//       int holesScoreWeight;
//       int outliersScoreWeight;
//       int otherScoreWeight;

//       int minHits;
//       int maxHoles;
//       int maxOutliers;
//       int maxUnused;
//       int maxSharedHits;

//       bool sharedHitsFlag;
//       std::size_t detectorId;

//       std::istringstream valueStream(value_str);
//       valueStream >> hitsScoreWeight ;
//       valueStream >> holesScoreWeight ;

//       valueStream >> outliersScoreWeight ;

//       valueStream >> otherScoreWeight ;

//       valueStream >> minHits ;

//       valueStream >> maxHoles ;

//       valueStream >> maxOutliers ;

//       valueStream >> maxUnused ;

//       valueStream >> maxSharedHits ;

//       valueStream >> sharedHitsFlag ;

//       valueStream >> detectorId ;

//       volumeMap[volumeId] = {
//         hitsScoreWeight, 
//         holesScoreWeight, 
//         outliersScoreWeight, 
//         otherScoreWeight,
//         minHits, 
//         maxHoles, 
//         maxOutliers, 
//         maxUnused, 
//         maxSharedHits, 
//         sharedHitsFlag, 
//         detectorId
//       };
//     }
//   }
//   return volumeMap;
// }



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

