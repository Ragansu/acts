// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/AmbiguityConfigJsonConverter.hpp"

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

#include <fstream>

namespace Acts {

std::pair<std::map<std::size_t, std::size_t>,
          std::map<std::size_t, AthenaAmbiguityResolution::DetectorConfig>>
AmbiguityConfigJsonConverter::fromJson(const std::string& configFile) const {
  std::ifstream file(configFile);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << configFile << std::endl;
    return {};
  }

  std::cout << "Reading configuration file: " << configFile << std::endl;

  nlohmann::json j;

  file >> j;

  std::map<std::size_t, AthenaAmbiguityResolution::DetectorConfig> detectorMap;
  std::map<std::size_t, std::size_t> volumeMap;

  for (auto& [key, value] : j.items()) {
    std::size_t detectorId = std::stoi(key);

    int hitsScoreWeight = value["hitsScoreWeight"];
    int holesScoreWeight = value["holesScoreWeight"];
    int outliersScoreWeight = value["outliersScoreWeight"];
    int otherScoreWeight = value["otherScoreWeight"];

    std::size_t minHits = value["minHits"];
    std::size_t maxHits = value["maxHits"];
    std::size_t maxHoles = value["maxHoles"];
    std::size_t maxDoubleHoles = value["maxDoubleHoles"];
    std::size_t maxOutliers = value["maxOutliers"];
    std::size_t maxSharedHits = value["maxSharedHits"];

    bool sharedHitsFlag = value["sharedHitsFlag"];

    std::vector<double> factorHits = value["factorHits"];
    std::vector<double> factorHoles = value["factorHoles"];

    auto detectorConfig = AthenaAmbiguityResolution::DetectorConfig(
        hitsScoreWeight, holesScoreWeight, outliersScoreWeight,
        otherScoreWeight, minHits, maxHits, maxDoubleHoles, maxHoles,
        maxOutliers, maxSharedHits, sharedHitsFlag, detectorId, factorHits,
        factorHoles);

    detectorMap[detectorId] = detectorConfig;

    std::vector<std::size_t> volumesIds = value["volumesIds"];
    for (auto volumeId : volumesIds) {
      volumeMap[volumeId] = detectorId;
    }
  }

  return std::make_pair(volumeMap, detectorMap);
}

}  // namespace Acts
