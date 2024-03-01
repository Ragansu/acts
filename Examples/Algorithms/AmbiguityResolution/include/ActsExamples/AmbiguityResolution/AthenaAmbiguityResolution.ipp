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
  std::vector<std::vector<std::tuple<std::size_t, std::size_t, Acts::ConstTrackStateType>>> 
    ActsExamples::AthenaAmbiguityResolution::computeInitialState(
      const ActsExamples::ConstTrackContainer& tracks,
      source_link_hash_t&& sourceLinkHash,
      source_link_equality_t&& sourceLinkEquality) const {
  auto measurementIndexMap =
      std::unordered_map<Acts::SourceLink, std::size_t, source_link_hash_t,
                         source_link_equality_t>(0, sourceLinkHash,
                                                 sourceLinkEquality);

  // std::vector<std::vector<std::size_t>> measurementsPerTrack;
  // std::vector<std::vector<std::size_t>> volumeIdsPerMeasurement;

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, Acts::ConstTrackStateType>>> measurementsPerTrack;


  int numberOfTracks = 0;
  for (const auto& track : tracks) {

    std::vector<std::tuple<std::size_t, std::size_t, Acts::ConstTrackStateType>> measurements_tuples;

    for (auto ts : track.trackStatesReversed()) {
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();

        const auto& geoID = ts.referenceSurface().geometryId();   

        // assign a new measurement index if the source link was not seen yet
        auto emplace = measurementIndexMap.try_emplace(
            sourceLink, measurementIndexMap.size());

        auto typeFlags = ts.typeFlags();
        
        measurements_tuples.push_back(std::make_tuple(emplace.first->second, geoID.volume(), typeFlags));
      }
    }

    measurementsPerTrack.push_back(std::move(measurements_tuples));

    ++numberOfTracks;
  }

  return measurementsPerTrack;
}
