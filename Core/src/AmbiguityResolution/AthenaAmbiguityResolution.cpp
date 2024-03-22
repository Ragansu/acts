// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/Measurement.hpp"

#include <stdexcept>

void Acts::AthenaAmbiguityResolution::VolumeConfig::setupScoreModifiers() {

  if (goodHits.size() != fakeHits.size()) {
      throw std::runtime_error("The number of good and fake hits must be the same");
      return;
  }
  else if(goodHits.size() != maxHits) {
    throw std::runtime_error("The number of goodHits = maxHits+1. Size of goodHits: " + std::to_string(goodHits.size()) + " maxHits: " + std::to_string(maxHits));
    return;
  } 
  for (std::size_t i=0; i< maxHits; ++i) m_factorHits.push_back(goodHits[i]/fakeHits[i]);

  if (goodHoles.size() != fakeHoles.size()){
    throw std::runtime_error("The number of good and fake holes must be the same");
    return;
  }
  else if(goodHoles.size() != maxHoles) {
    throw std::runtime_error("The number of goodHoles = maxHoles+1. Size of goodHoles: " + std::to_string(goodHoles.size()) + " maxHoles: " + std::to_string(maxHoles));
    return;
  }
  for (std::size_t i=0; i<maxHoles; ++i) m_factorHoles.push_back(goodHoles[i]/fakeHoles[i]);
}

std::vector<std::size_t> Acts::AthenaAmbiguityResolution::getCleanedOutTracks(
    std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack) const {
    std::vector<std::size_t> cleanTracks;

  ACTS_INFO("Cleaning tracks");
  
  std::size_t numberOfTracks = measurementsPerTrack.size();
  std::map<std::size_t,VolumeConfig> detectorMap;
  std::vector<std::size_t> detectorList;
  std::cout << "Number of tracks: " << numberOfTracks << std::endl;

  boost::container::flat_map<std::size_t,
                            boost::container::flat_set<std::size_t>>
    tracksPerMeasurement;

  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      auto iVolume = std::get<1>(measurements_tuples);
      auto detector_it = m_cfg.volumeMap.find(iVolume);
      if(detector_it != m_cfg.volumeMap.end()){
        auto detector = detector_it->second;

        auto volume_it = detectorMap.find(detector.detectorId);
        if(volume_it == detectorMap.end()){
          detectorMap[detector.detectorId] = detector;
          detectorList.push_back(detector.detectorId);
        }      
        tracksPerMeasurement[iMeasurement].insert(iTrack);
      }
    }
  }

  enum TsosTypes {
    // A measurement not yet used in any other track
    UnusedHit   = 1,
    // A measurement shared with another track
    SharedHit   = 2,
    // A hit that needs to be removed from the track
    RejectedHit = 3,
    // an outlier, to be copied in case
    Outlier     = 4,
    // other TSOS types to be copied in case
    OtherTsos   = 5
  };

 
  std::vector<std::vector<std::size_t>>  newMeasurements;
  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {

    int track_score = trackScore[iTrack];
    std::cout << "Track score: " << track_score << std::endl;


    int numUnused         = 0;
    int numShared         = 0;

    // int numWeightedShared = 0;

    auto counterMap = counterMaps[iTrack];

    bool TrkCouldBeAccepted = true;

    std::cout << "Number of Detectors: " << detectorList.size() << std::endl;

    // for tracks with shared hits, we need to check and remove bad hits

    std::vector<int> tsosTypes(measurementsPerTrack[iTrack].size(), OtherTsos);
    int index = 0;


    // std::cout << " computing tsos types for track " << iTrack << std::endl;
    // std::cout << " Number of measurements: " << measurementsPerTrack[iTrack].size() << std::endl;

    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {

      tsosTypes[index] = OtherTsos;
    
      auto iMeasurement = std::get<0>(measurements_tuples);
      auto iVolume = std::get<1>(measurements_tuples);
      auto isoutliner = std::get<2>(measurements_tuples);

      auto detector_it = m_cfg.volumeMap.find(iVolume);
      if(detector_it == m_cfg.volumeMap.end()){
        index++;
        continue;
      }

      auto detector = detector_it->second;

      if (isoutliner) {
        ACTS_INFO ("Measurement is outlier on a fitter track, copy it over");
        tsosTypes[index] = Outlier;
        index++;
        continue;
      }  
      if (tracksPerMeasurement[iMeasurement].size() == 1) {
        ACTS_VERBOSE ("Measurement is not shared, copy it over");

        tsosTypes[index] = UnusedHit;
        (counterMap[detector.detectorId].nUnused)++;
        numUnused++;

        index++;
        continue;
      }        
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
          ACTS_VERBOSE ("Measurement is shared, copy it over");

        if (detector.sharedHitsFlag == true) {
          ACTS_VERBOSE ("Measurement is shared, Reject it");
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
    std::size_t measurment = 0;
    std::size_t cntIns = 0;

    for( std::size_t i = 0; i < tsosTypes.size(); i++){

      auto measurement_tuples  = measurementsPerTrack[iTrack][i];
      measurment = std::get<0>(measurement_tuples);
      
      std::cout << "Tsos type " << i << " is " << tsosTypes[i] << std::endl;
      if (tsosTypes[i] == RejectedHit){
        std::cout << "Droping rejected hit" << std::endl;
      }
      else if (tsosTypes[i] != SharedHit){
        std::cout << "Good TSOS, copy hit" << std::endl;
        newMeasurementsPerTrack.push_back(measurment);
      }
      else if (cntIns >= m_cfg.maxShared){
          std::cout << "Too many shared hit, drop it" << std::endl;
      }
      else {
        std::cout << "Try to recover shared hit " << std::endl;
        if (tracksPerMeasurement[measurment].size() < m_cfg.maxSharedTracksPerMeasurement &&
            track_score > m_cfg.minScoreSharedTracks) {
          
          std::cout << "Accepted hit shared with " << tracksPerMeasurement[measurment].size() << " tracks" <<std::endl;
          newMeasurementsPerTrack.push_back(measurment);
          cntIns++;
        }
        else {
          std::cout << "Rejected hit shared with " << tracksPerMeasurement[measurment].size() << " tracks" <<std::endl;
        }
      }
    }

    if (newMeasurementsPerTrack.size() < 3) {
      TrkCouldBeAccepted = false;
      std::cout << "Track " << iTrack << " could not be accepted - not enought hits" << std::endl;
      continue;
    }
    
    for(long unsigned int i = 0; i< detectorList.size(); i++){

      auto volume_it = detectorMap.find(detectorList[i]);
      auto detector = volume_it->second;
      if (counterMap[detector.detectorId].nSharedHits > detector.maxSharedHits){
        TrkCouldBeAccepted = false;
        break;
      }
    }

    if (TrkCouldBeAccepted){
      cleanTracks.push_back(iTrack);
      newMeasurements.push_back(newMeasurementsPerTrack);
      continue;
    }
  }

  numberOfTracks = cleanTracks.size();


  boost::container::flat_map<std::size_t,boost::container::flat_set<std::size_t>> newTracksPerMeasurement;

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
    std::cout << "Number of shared measurements: " << sharedMeasurementsPerTrack << std::endl;
    std::cout << "Number of measurements: " << newMeasurements[track_id].size() << std::endl;
    std::cout << "Score of the track: " << trackScore[cleanTracks[track_id]] << std::endl;
    std::cout << "------------------------------------ " << std::endl;
  }

  return cleanTracks;
}