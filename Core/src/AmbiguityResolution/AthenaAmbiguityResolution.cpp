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
    std::vector<int> trackScore,
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
  std::cout << "Number of tracks: " << numberOfTracks << std::endl;

  boost::container::flat_map<std::size_t,
                             boost::container::flat_set<std::size_t>>
      tracksPerMeasurement;

  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    if (trackScore[iTrack] <= 0) {
      measurementsPerTrack[iTrack].clear();
      continue;
    }
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      tracksPerMeasurement[iMeasurement].insert(iTrack);
      std::cout << "Track " << iTrack << " has measurement " << iMeasurement
                << std::endl;
    }
  }

  enum TsosTypes {
    // A measurement not yet used in any other track
    UnusedHit = 1,
    // A measurement shared with another track
    SharedHit = 2,
    // A hit that needs to be removed from the track
    RejectedHit = 3,
    // an outlier, to be copied in case
    Outlier = 4,
    // other TSOS types to be copied in case
    OtherTsos = 5
  };

  std::vector<std::vector<std::size_t>> newMeasurements;

  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    int track_score = trackScore[iTrack];
    std::cout << "Track score: " << track_score << std::endl;

    if (track_score <= 0) {
      std::cout << "Track " << iTrack << " could not be accepted - low score"
                << std::endl;
      continue;
    }

    int numUnused = 0;
    int numShared = 0;

    // int numWeightedShared = 0;

    auto counterMap = counterMaps[iTrack];

    bool TrkCouldBeAccepted = true;

    std::cout << "Number of Detectors: " << m_cfg.detectorMap.size()
              << std::endl;

    // for tracks with shared hits, we need to check and remove bad hits

    std::vector<int> tsosTypes(measurementsPerTrack[iTrack].size(), OtherTsos);
    int index = 0;

    // std::cout << " computing tsos types for track " << iTrack << std::endl;
    // std::cout << " Number of measurements: " <<
    // measurementsPerTrack[iTrack].size() << std::endl;

    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      tsosTypes[index] = OtherTsos;

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
        tsosTypes[index] = Outlier;
        index++;
        continue;
      }
      if (tracksPerMeasurement[iMeasurement].size() == 1) {
        ACTS_VERBOSE("Measurement is not shared, copy it over");

        tsosTypes[index] = UnusedHit;
        (counterMap[detectorId].nUnused)++;
        numUnused++;

        index++;
        continue;
      }
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
        ACTS_VERBOSE("Measurement is shared, copy it over");

        if (detector.sharedHitsFlag == true) {
          ACTS_VERBOSE("Measurement is shared, Reject it");
          tsosTypes[index] = RejectedHit;
          index++;
          continue;
        }

        tsosTypes[index] = SharedHit;
        numShared++;

        // Yet to be implemented
        // numWeightedShared += (isPixel ? 2 : 1);

        index++;
        continue;
      }
    }

    std::vector<std::size_t> newMeasurementsPerTrack;
    std::size_t measurement = 0;
    std::size_t cntIns = 0;

    for (std::size_t i = 0; i < tsosTypes.size(); i++) {
      auto measurement_tuples = measurementsPerTrack[iTrack][i];
      measurement = std::get<0>(measurement_tuples);

      // std::cout << "Tsos type " << i << " is " << tsosTypes[i] << std::endl;
      if (tsosTypes[i] == RejectedHit) {
        // std::cout << "Dropping rejected hit" << std::endl;
      } else if (tsosTypes[i] != SharedHit) {
        // std::cout << "Good TSOS, copy hit" << std::endl;
        newMeasurementsPerTrack.push_back(measurement);
      } else if (cntIns >= m_cfg.maxShared) {
        // std::cout << "Too many shared hit, drop it" << std::endl;
      } else {
        // std::cout << "Try to recover shared hit " << std::endl;
        if (tracksPerMeasurement[measurement].size() <
                m_cfg.maxSharedTracksPerMeasurement &&
            track_score > m_cfg.minScoreSharedTracks) {
          std::cout << "Accepted hit shared with "
                    << tracksPerMeasurement[measurement].size() << " tracks"
                    << std::endl;
          newMeasurementsPerTrack.push_back(measurement);
          cntIns++;
        } else {
          std::cout << "Rejected hit shared with "
                    << tracksPerMeasurement[measurement].size() << " tracks"
                    << std::endl;
        }
      }
    }

    if (newMeasurementsPerTrack.size() < 3) {
      TrkCouldBeAccepted = false;
      std::cout << "Track " << iTrack
                << " could not be accepted - not enough hits" << std::endl;
      std::cout << "Number of hits: " << measurementsPerTrack[iTrack].size()
                << std::endl;
      std::cout << "Number of good hits: " << newMeasurementsPerTrack.size()
                << std::endl;
      std::cout << "-------------------------------------" << std::endl;
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
      std::cout << "Track " << iTrack << " is accepted" << std::endl;
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
    std::cout << "Cleaned track " << cleanTracks[track_id] << std::endl;
    std::size_t sharedMeasurementsPerTrack = 0;
    for (auto iMeasurement : newMeasurements[track_id]) {
      if (newTracksPerMeasurement[iMeasurement].size() > 1) {
        ++sharedMeasurementsPerTrack;
      }
    }
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Track ID: " << cleanTracks[track_id] << std::endl;
    std::cout << "Number of shared measurements: " << sharedMeasurementsPerTrack
              << std::endl;
    std::cout << "Number of measurements: " << newMeasurements[track_id].size()
              << std::endl;
    std::cout << "Score of the track: " << trackScore[cleanTracks[track_id]]
              << std::endl;
    std::cout << "------------------------------------ " << std::endl;
  }

  return cleanTracks;
}
