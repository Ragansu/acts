// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolution.hpp"

#include <unordered_map>


template <typename source_link_hash_t,
          typename source_link_equality_t>
std::vector<std::vector<std::size_t>> ActsExamples::AthenaAmbiguityResolution::computeInitialState(
    const ActsExamples::ConstTrackContainer& tracks,
    source_link_hash_t&& sourceLinkHash,
    source_link_equality_t&& sourceLinkEquality) const {
  auto measurementIndexMap =
      std::unordered_map<Acts::SourceLink, std::size_t, source_link_hash_t,
                         source_link_equality_t>(0, sourceLinkHash,
                                                 sourceLinkEquality);

  // Iterate through all input tracks, collect their properties like measurement
  // count and chi2 and fill the measurement map in order to relate tracks to
  // each other if they have shared hits.

  // std::vector<std::size_t> sharedMeasurementsPerTrack;
  std::vector<std::vector<std::size_t>> measurementsPerTrack;
  // std::vector<std::set<std::size_t>> tracksPerMeasurement;
  int numberOfTracks = 0;
  for (const auto& track : tracks) {
    // Kick out tracks that do not fulfill our initial requirements

    std::vector<std::size_t> measurements;

    for (auto ts : track.trackStatesReversed()) {
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();
        // assign a new measurement index if the source link was not seen yet
        auto emplace = measurementIndexMap.try_emplace(
            sourceLink, measurementIndexMap.size());
        measurements.push_back(emplace.first->second);
      }
    }

    measurementsPerTrack.push_back(std::move(measurements));

    ++numberOfTracks;
  }

  return measurementsPerTrack;
}
