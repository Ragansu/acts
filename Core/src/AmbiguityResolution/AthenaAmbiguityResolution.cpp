// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <stdexcept>

std::vector<std::size_t> Acts::AthenaAmbiguityResolution::getCleanedOutTracks(
    std::vector<double> trackScore,
    std::vector<std::map<std::size_t, Counter>>& counterMaps,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
        measurementsPerTrack) const {
  std::vector<std::size_t> cleanTracks;

  ACTS_INFO("Cleaning tracks");

  if (trackScore.size() != measurementsPerTrack.size()) {
    throw std::invalid_argument(
        "Track score and measurementsPerTrack size mismatch");
  }

  std::size_t numberOfTracks = measurementsPerTrack.size();
  ACTS_DEBUG("Number of tracks: " << numberOfTracks);

  boost::container::flat_map<std::size_t,
                             boost::container::flat_set<std::size_t>>
      tracksPerMeasurement;

  // Removes bad tracks and counts computes the vector of tracks per
  // measurement.
  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    if (trackScore[iTrack] <= 0) {
      measurementsPerTrack[iTrack].clear();
      continue;
    }
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      tracksPerMeasurement[iMeasurement].insert(iTrack);
    }
  }

  enum TrackStateTypes {
    // A measurement not yet used in any other track
    UnusedHit = 1,
    // A measurement shared with another track
    SharedHit = 2,
    // A hit that needs to be removed from the track
    RejectedHit = 3,
    // an outlier, to be copied in case
    Outlier = 4,
    // other trackstate types to be copied in case
    OtherTrackStateType = 5
  };

  std::vector<std::vector<std::size_t>> newMeasurements;

  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    double track_score = trackScore[iTrack];
    ACTS_DEBUG("Track score: " << track_score);

    if (track_score <= 0) {
      ACTS_DEBUG("Track " << iTrack << " could not be accepted - low score");
      continue;
    }

    auto counterMap = counterMaps[iTrack];

    bool TrkCouldBeAccepted = true;

    // for tracks with shared hits, we need to check and remove bad hits

    std::vector<int> trackStateTypes(measurementsPerTrack[iTrack].size(),
                                     OtherTrackStateType);
    int index = 0;

    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      auto iVolume = std::get<1>(measurements_tuples);
      auto isoutliner = std::get<2>(measurements_tuples);

      auto volume_it = m_cfg.volumeMap.find(iVolume);

      if (volume_it == m_cfg.volumeMap.end()) {
        index++;
        continue;
      }

      auto detectorId = volume_it->second;
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;

      if (isoutliner) {
        ACTS_INFO("Measurement is outlier on a fitter track, copy it over");
        trackStateTypes[index] = Outlier;
        index++;
        continue;
      }
      if (tracksPerMeasurement[iMeasurement].size() == 1) {
        ACTS_VERBOSE("Measurement is not shared, copy it over");

        trackStateTypes[index] = UnusedHit;

        index++;
        continue;
      }
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
        ACTS_VERBOSE("Measurement is shared, copy it over");

        if (detector.sharedHitsFlag == true) {
          ACTS_VERBOSE("Measurement is shared, Reject it");
          trackStateTypes[index] = RejectedHit;
          index++;
          continue;
        }

        trackStateTypes[index] = SharedHit;

        index++;
        continue;
      }
    }

    std::vector<std::size_t> newMeasurementsPerTrack;
    std::size_t measurement = 0;
    std::size_t cntIns = 0;

    for (std::size_t i = 0; i < trackStateTypes.size(); i++) {
      auto measurement_tuples = measurementsPerTrack[iTrack][i];
      measurement = std::get<0>(measurement_tuples);

      if (trackStateTypes[i] == RejectedHit) {
        ACTS_DEBUG("Dropping rejected hit");
      } else if (trackStateTypes[i] != SharedHit) {
        ACTS_DEBUG("Good TSOS, copy hit");
        newMeasurementsPerTrack.push_back(measurement);
      } else if (cntIns >= m_cfg.maxShared) {
        ACTS_DEBUG("Too many shared hit, drop it");
      } else {
        ACTS_DEBUG("Try to recover shared hit ");
        if (tracksPerMeasurement[measurement].size() <
                m_cfg.maxSharedTracksPerMeasurement &&
            track_score > m_cfg.minScoreSharedTracks) {
          ACTS_DEBUG("Accepted hit shared with "
                     << tracksPerMeasurement[measurement].size() << " tracks");
          newMeasurementsPerTrack.push_back(measurement);
          cntIns++;
        } else {
          ACTS_DEBUG("Rejected hit shared with "
                     << tracksPerMeasurement[measurement].size() << " tracks");
        }
      }
    }

    if (newMeasurementsPerTrack.size() < 3) {
      TrkCouldBeAccepted = false;
      ACTS_DEBUG(std::endl
                 << "Track " << iTrack
                 << " could not be accepted - not enough hits");
      ACTS_DEBUG("Number of hits: " << measurementsPerTrack[iTrack].size());
      ACTS_DEBUG("Number of good hits: " << newMeasurementsPerTrack.size());
      continue;
    }

    for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
         detectorId++) {
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;
      if (counterMap[detectorId].nSharedHits > detector.maxSharedHits) {
        TrkCouldBeAccepted = false;
        break;
      }
    }

    if (TrkCouldBeAccepted) {
      cleanTracks.push_back(iTrack);
      newMeasurements.push_back(newMeasurementsPerTrack);
      ACTS_VERBOSE("Track " << iTrack << " is accepted");
      continue;
    }
  }

  numberOfTracks = cleanTracks.size();

  boost::container::flat_map<std::size_t,
                             boost::container::flat_set<std::size_t>>
      newTracksPerMeasurement;

  for (std::size_t track_id = 0; track_id < numberOfTracks; ++track_id) {
    for (auto iMeasurement : newMeasurements[track_id]) {
      newTracksPerMeasurement[iMeasurement].insert(track_id);
    }
  }

  for (std::size_t track_id = 0; track_id < numberOfTracks; ++track_id) {
    std::size_t sharedMeasurementsPerTrack = 0;
    for (auto iMeasurement : newMeasurements[track_id]) {
      if (newTracksPerMeasurement[iMeasurement].size() > 1) {
        ++sharedMeasurementsPerTrack;
      }
    }
    ACTS_VERBOSE(std::endl << "Track ID: " << cleanTracks[track_id]);
    ACTS_VERBOSE(
        "Number of shared measurements: " << sharedMeasurementsPerTrack);
    ACTS_VERBOSE("Number of measurements: " << newMeasurements[track_id].size())
    ACTS_VERBOSE("Score of the track: " << trackScore[cleanTracks[track_id]]);
  }

  return cleanTracks;
}
