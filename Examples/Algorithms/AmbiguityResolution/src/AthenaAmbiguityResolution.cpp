// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

ActsExamples::AthenaAmbiguityResolution::AthenaAmbiguityResolution(
    std::string name, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm(name, lvl) {}


ActsExamples::ConstTrackContainer
ActsExamples::AthenaAmbiguityResolution::prepareOutputTrack(
    const ActsExamples::ConstTrackContainer& tracks,
    std::vector<std::size_t>& goodTracks) const {
  std::shared_ptr<Acts::ConstVectorMultiTrajectory> trackStateContainer =
      tracks.trackStateContainerHolder();
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  trackContainer->reserve(goodTracks.size());
  // Temporary empty track state container: we don't change the original one,
  // but we need one for filtering
  auto tempTrackStateContainer =
      std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer solvedTracks{trackContainer, tempTrackStateContainer};
  solvedTracks.ensureDynamicColumns(tracks);

  for (auto&& iTrack : goodTracks) {
    auto destProxy = solvedTracks.getTrack(solvedTracks.addTrack());
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFrom(srcProxy, false);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      trackStateContainer};
  return outputTracks;
}

std::vector<int> ActsExamples::AthenaAmbiguityResolution::simpleScore(
 const ActsExamples::ConstTrackContainer& tracks, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {

  std::vector<int> trackScore;
  int iTrack = 0;  


  // Loop over all the trajectories in the events
  for (auto track : tracks){
    auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
      tracks.trackStateContainer(), track.tipIndex());
    int score = 100;
    auto counterMap = m_counterMap;
    
    // if (track.chi2() > 0 && track.nDoF() > 0) {
    //   score+= std::log10(1.0-(track.chi2()/track.nDoF())); // place holder
    // }

    // TODO: add scored based on chi2 and ndof

    // detector score is determined by the number of hits/hole/outliers * hit/hole/outlier score
    // here so instead of calculating nhits/nholes/noutliers per volume, 
    // we just loop over all volumes and add the score.

    for (long unsigned int i = 0; i < trajState.measurementVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.measurementVolume[i]);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;
        score+=detector.hitsScoreWeight;
        counterMap[detector.detectorId].nhits++;
      }
    }

    for (long unsigned int i = 0; i < trajState.holeVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.holeVolume[i]);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;
        score+=detector.holesScoreWeight;
        counterMap[detector.detectorId].nholes++;
      }
      else{
        ACTS_INFO("Detector not found");
        ACTS_INFO("Detector ID: " << trajState.holeVolume[i]);
      } 
    }
    for (long unsigned int i = 0; i < trajState.outlierVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.outlierVolume[i]);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;
        score+=detector.outliersScoreWeight;
        counterMap[detector.detectorId].noutliers++;
      }
    }
      // TODO: add scored based on eta and phi

    trackScore.push_back(score);
    // ACTS_INFO("Track " << iTrack << " score: " << score);

    // measurementPerTrack.push_back(track.measurements();

    counterMaps.push_back(counterMap);
    iTrack++;


  } // end of loop over tracks
    
  return trackScore;
}

std::size_t sourceLinkHash(const Acts::SourceLink& a) {
  return static_cast<std::size_t>(
      a.get<ActsExamples::IndexSourceLink>().index());
}

bool sourceLinkEquality(const Acts::SourceLink& a, const Acts::SourceLink& b) {
  return a.get<ActsExamples::IndexSourceLink>().index() ==
         b.get<ActsExamples::IndexSourceLink>().index();
}

// place holder for goodTracks algorithm
std::vector<std::size_t> 
ActsExamples::AthenaAmbiguityResolution::solveAmbiguity(
    const ActsExamples::ConstTrackContainer& tracks ,std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {

  std::vector<std::size_t> cleanTracks = getCleanedOutTracks(tracks, trackScore, counterMaps);

  ACTS_INFO("Number of tracks: " << tracks.size());

  ACTS_INFO("Number of clean tracks: " << cleanTracks.size());

  std::vector<std::size_t> goodTracks;

  ACTS_INFO("Min score: " << m_minScore);


  for(long unsigned int i=0; i<cleanTracks.size(); ++i){
    // ACTS_INFO("Track " << i << " score: " << trackScore[cleanTracks[i]]);
    if (trackScore[cleanTracks[i]] > m_minScore){
      goodTracks.push_back(cleanTracks[i]);
    }
  }
  ACTS_INFO("Number of good tracks: " << goodTracks.size());
  return goodTracks;
}


std::vector<std::size_t> ActsExamples::AthenaAmbiguityResolution::getCleanedOutTracks(
    const ActsExamples::ConstTrackContainer& tracks ,std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {
  std::vector<std::size_t> cleanTracks;

  // Loop over all detectors

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, Acts::ConstTrackStateType>>>  measurementsPerTrack = computeInitialState(tracks, &sourceLinkHash, &sourceLinkEquality);

  boost::container::flat_map<std::size_t,boost::container::flat_set<std::size_t>> tracksPerMeasurement;
  std::size_t numberOfTracks = tracks.size();



  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      tracksPerMeasurement[iMeasurement].insert(iTrack);
    }
  }

  std::vector<std::size_t> sharedMeasurementsPerTrack = std::vector<std::size_t>(tracks.size(), 0);
  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
        ++sharedMeasurementsPerTrack[iTrack];
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


  int iTrack = 0;
  for (const auto& track : tracks) {

    int numUnused         = 0;
    int numShared         = 0;

    // int numWeightedShared = 0;

    auto counterMap = counterMaps[iTrack];

    auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(tracks.trackStateContainer(), track.tipIndex());
    bool TrkCouldBeAccepted = true;
    auto volumeList = trajState.measurementVolume;
    auto volumeIterator = std::unique(volumeList.begin(), volumeList.end());

    volumeList.resize(std::distance(volumeList.begin(), volumeIterator));
 
    for(long unsigned int i = 0; i< volumeList.size(); ++i){
      auto detector_it = m_volumeMap.find(volumeList[i]);
      if(detector_it == m_volumeMap.end()){
        continue;
      }
      auto detector = detector_it->second;
      
      ACTS_DEBUG ("---> Found summary information");
      ACTS_DEBUG ("---> Detector ID: " << detector.detectorId);
      ACTS_DEBUG ("---> Number of hits: " << counterMap[detector.detectorId].nhits);
      ACTS_DEBUG ("---> Number of holes: " << counterMap[detector.detectorId].nholes);
      ACTS_DEBUG ("---> Number of outliers: " << counterMap[detector.detectorId].noutliers);

      if (counterMap[detector.detectorId].nhits < detector.minHits){
        TrkCouldBeAccepted = false;
        break;
      }

      if (counterMap[detector.detectorId].nholes > detector.maxHoles){
        TrkCouldBeAccepted = false;
        break;
      }

      if (counterMap[detector.detectorId].noutliers > detector.maxOutliers){
        TrkCouldBeAccepted = false;
        break;
      }
      
    }

    if (!TrkCouldBeAccepted){
      iTrack++;
      continue;
    }

    // for tracks with shared hits, we need to check and remove bad hits

    std::vector<int> tsosTypes(trajState.measurementVolume.size(), 0);
    int index = 0;
    bool firstisshared = true;

    std::size_t* lastrot       = nullptr;
    std::size_t* lastbutonerot = nullptr;
    int          lastrotindex  = 0;

    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {

      tsosTypes[index] = OtherTsos;
    
      auto iMeasurement = std::get<0>(measurements_tuples);
      auto iVolume = std::get<1>(measurements_tuples);
      auto iTypeFlags = std::get<2>(measurements_tuples);

      auto detector_it = m_volumeMap.find(iVolume);
      if(detector_it == m_volumeMap.end()){
        continue;
      }

      auto detector = detector_it->second;

      if (iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
        ACTS_VERBOSE ("Measurement is outlier on a fitter track, copy it over");
        tsosTypes[index] = Outlier;
        continue;
      }  
          
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
          ACTS_VERBOSE ("Measurement is shared, copy it over");
        
        if (!iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
          ACTS_VERBOSE ("Measurement is outlier on a shared track, copy it over");

          tsosTypes[index] = UnusedHit;
          (counterMap[detector.detectorId].nUnused)++;
          numUnused++;
          if (numShared == 0) firstisshared = false;

          lastbutonerot = lastrot;
          lastrot       = &iMeasurement;
          lastrotindex  = index;

          continue;
        }

        if (lastrot != nullptr && lastbutonerot != nullptr) {

          std::cout << lastrotindex << std::endl;
          std::cout << firstisshared << std::endl;

        }


        tsosTypes[index] = SharedHit;
        (counterMap[detector.detectorId].nSharedHits)++;
        numShared++;

        // Yet to be implemented
        // numWeightedShared += (isPixel ? 2 : 1);

        lastbutonerot = lastrot;
        lastrot       = &iMeasurement;
        lastrotindex  = index;    

        continue;      
      }

      ACTS_VERBOSE ("Measurement is not shared, Reject it");
      tsosTypes[index] = RejectedHit;
      TrkCouldBeAccepted = false;
      index++;
    }

    if (sharedMeasurementsPerTrack[iTrack] > 0 && trackScore[iTrack] < m_minScoreSharedTracks) {
      TrkCouldBeAccepted = false;
      iTrack++;
      continue;
    }
    // temporary solution for access detector information
    for ( std::size_t i = 0; i< volumeList.size(); ++i){
      auto detector_it = m_volumeMap.find(volumeList[i]);
      if(detector_it == m_volumeMap.end()){
        continue;
      }
      auto detector = detector_it->second;
      if (counterMap[detector.detectorId].nSharedHits > detector.maxSharedHits){

        TrkCouldBeAccepted = false;
        iTrack++;
      }
    }


    if (TrkCouldBeAccepted){
      cleanTracks.push_back(iTrack);
      // ACTS_INFO("Track " << iTrack << " is clean");
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Track " << iTrack << " has " << sharedMeasurementsPerTrack[iTrack] << " shared measurements" << std::endl;
      std::cout << "Track " << iTrack << " has " << measurementsPerTrack[iTrack].size() << " measurements" << std::endl;
      std::cout << "Track " << iTrack << " has " << trackScore[iTrack] << " score" << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    iTrack++;

  }
  return cleanTracks;
}


