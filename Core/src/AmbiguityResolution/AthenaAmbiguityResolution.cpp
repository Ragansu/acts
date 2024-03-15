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



std::vector<std::size_t> Acts::AthenaAmbiguityResolution::getCleanedOutTracks(
    std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack) const {
    std::vector<std::size_t> cleanTracks;

  std::cout << "Cleaning tracks" << std::endl;
  
  std::size_t numberOfTracks = measurementsPerTrack.size();
  std::vector<std::size_t> volumeList; 

  std::cout << "Number of tracks: " << numberOfTracks << std::endl;

  boost::container::flat_map<std::size_t,
                            boost::container::flat_set<std::size_t>>
    tracksPerMeasurement;

  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      auto iVolume = std::get<1>(measurements_tuples);
      volumeList.push_back(iVolume);
      tracksPerMeasurement[iMeasurement].insert(iTrack);
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


    int numUnused         = 0;
    int numShared         = 0;

    // int numWeightedShared = 0;

    auto counterMap = counterMaps[iTrack];

    bool TrkCouldBeAccepted = true;

    // auto volumeList = trajState.measurementVolume;
    auto volumeIterator = std::unique(volumeList.begin(), volumeList.end());

    volumeList.resize(std::distance(volumeList.begin(), volumeIterator));

    std::cout << "Number of volumes: " << volumeList.size() << std::endl;
 
    for(long unsigned int i = 0; i< volumeList.size(); ++i){
      auto detector_it = m_cfg.volumeMap.find(volumeList[i]);
      if(detector_it == m_cfg.volumeMap.end()){
        continue;
      }
      auto detector = detector_it->second;
      
      ACTS_DEBUG ("---> Found summary information");
      ACTS_DEBUG ("---> Detector ID: " << detector.detectorId);
      ACTS_DEBUG ("---> Number of hits: " << counterMap[detector.detectorId].nHits);
      ACTS_DEBUG ("---> Number of holes: " << counterMap[detector.detectorId].nHoles);
      ACTS_DEBUG ("---> Number of outliers: " << counterMap[detector.detectorId].nOutliers);

      if (counterMap[detector.detectorId].nHits < detector.minHits){
        TrkCouldBeAccepted = false;
      }

      if (counterMap[detector.detectorId].nHoles > detector.maxHoles){
        TrkCouldBeAccepted = false;
      }

      if (counterMap[detector.detectorId].nOutliers > detector.maxOutliers){
        TrkCouldBeAccepted = false;
      }
      
    }

    if (!TrkCouldBeAccepted){
      std::cout << "Track " << iTrack << " could not be accepted" << std::endl;
      continue;
    }

    // for tracks with shared hits, we need to check and remove bad hits

    std::vector<int> tsosTypes(measurementsPerTrack[iTrack].size(), OtherTsos);
    int index = 0;


    // std::cout << " computing tsos types for track " << iTrack << std::endl;
    // std::cout << " Number of measurements: " << measurementsPerTrack[iTrack].size() << std::endl;

    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {

      // std::cout << " computing tsos for measurement " << index << std::endl;

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
        (counterMap[detector.detectorId].nSharedHits)++;
        numShared++;

        // Yet to be implemented
        // numWeightedShared += (isPixel ? 2 : 1);

         index++;
        continue;      
      }

      ACTS_VERBOSE ("Measurement is not shared, Reject it");
      tsosTypes[index] = RejectedHit;
      index++;
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
      else if (cntIns >= m_maxShared){
          std::cout << "Too many shared hit, drop it" << std::endl;
      }
      else {
        std::cout << "Try to recover shared hit " << std::endl;
        if (tracksPerMeasurement[measurment].size() < m_maxSharedTracksPerMeasurement &&
            track_score > m_minScoreSharedTracks) {
          
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
    
    for ( std::size_t i = 0; i< volumeList.size(); ++i){
      auto detector_it = m_cfg.volumeMap.find(volumeList[i]);
      if(detector_it == m_cfg.volumeMap.end()){
        continue;
      }
      auto detector = detector_it->second;
      if (counterMap[detector.detectorId].nSharedHits > detector.maxSharedHits){
        TrkCouldBeAccepted = false;
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
    std::cout << "Score of the track: " << trackScore[track_id] << std::endl;
    std::cout << "------------------------------------ " << std::endl;
  }

  return cleanTracks;
}