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

std::size_t sourceLinkHash(const Acts::SourceLink& a) {
  return static_cast<std::size_t>(
      a.get<ActsExamples::IndexSourceLink>().index());
}

bool sourceLinkEquality(const Acts::SourceLink& a, const Acts::SourceLink& b) {
  return a.get<ActsExamples::IndexSourceLink>().index() ==
         b.get<ActsExamples::IndexSourceLink>().index();
}

std::vector<int> ActsExamples::AthenaAmbiguityResolution::simpleScore(
 const ActsExamples::ConstTrackContainer& tracks, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {

  std::vector<int> trackScore;
  int iTrack = 0;  

  // Loop over all the trajectories in the events
  for (auto track : tracks){

    int score = 100;
    auto counterMap = m_counterMap;
    
    if (track.chi2() > 0 && track.nDoF() > 0) {
      score+= (1.0-std::log10(track.chi2()/track.nDoF())); // place holder
    }

    // TODO: add scored based on chi2 and ndof


    // detector score is determined by the number of hits/hole/outliers * hit/hole/outlier score
    // here so instead of calculating nHits/nHoles/nOutliers per volume, 
    // we just loop over all volumes and add the score.


    for (auto ts : track.trackStatesReversed()) {
      auto iVolume = ts.referenceSurface().geometryId().volume();
      auto iTypeFlags = ts.typeFlags();  

      auto detector_it = m_volumeMap.find(iVolume);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;

        if (iTypeFlags.test(Acts::TrackStateFlag::MeasurementFlag)){
          if (iTypeFlags.test(Acts::TrackStateFlag::SharedHitFlag)){
            counterMap[detector.detectorId].nSharedHits++;
          }
          score+=detector.hitsScoreWeight;
          counterMap[detector.detectorId].nHits++;
        }
        else if (iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
          score+=detector.outliersScoreWeight;
          counterMap[detector.detectorId].nOutliers++;
        }
        else if (iTypeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
          score+=detector.holesScoreWeight;
          counterMap[detector.detectorId].nHoles++;
        }
      }
    }


    trackScore.push_back(score);
    ACTS_INFO("Track " << iTrack << " score: " << score);


    counterMaps.push_back(counterMap);
    iTrack++;


  } // end of loop over tracks
    
  return trackScore;
}


// place holder for goodTracks algorithm
std::vector<std::size_t> 
ActsExamples::AthenaAmbiguityResolution::solveAmbiguity(
    const ActsExamples::ConstTrackContainer& tracks ,std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {

  std::cout << "Solving ambiguity" << std::endl;

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

  std::cout << "Cleaning tracks" << std::endl;

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, Acts::ConstTrackStateType>>>  measurementsPerTrack = computeInitialState(tracks, &sourceLinkHash, &sourceLinkEquality);

  boost::container::flat_map<std::size_t,boost::container::flat_set<std::size_t>> tracksPerMeasurement;
  
  std::cout << "Number of measurements: " << measurementsPerTrack.size() << std::endl;
  
  std::size_t numberOfTracks = tracks.size();


  std::cout << "Number of tracks: " << numberOfTracks << std::endl;


  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurements_tuples);
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

    std::cout << "Number of volumes: " << volumeList.size() << std::endl;
 
    for(long unsigned int i = 0; i< volumeList.size(); ++i){
      auto detector_it = m_volumeMap.find(volumeList[i]);
      if(detector_it == m_volumeMap.end()){
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

    // std::cout << " computing tsos types for track " << iTrack << std::endl;
    // std::cout << " Number of measurements: " << measurementsPerTrack[iTrack].size() << std::endl;

    for (auto measurements_tuples : measurementsPerTrack[iTrack]) {

      // std::cout << " computing tsos for measurement " << index << std::endl;

      tsosTypes[index] = OtherTsos;
    
      auto iMeasurement = std::get<0>(measurements_tuples);
      auto iVolume = std::get<1>(measurements_tuples);
      auto iTypeFlags = std::get<2>(measurements_tuples);

      auto detector_it = m_volumeMap.find(iVolume);
      if(detector_it == m_volumeMap.end()){
        index++;
        continue;
      }

      auto detector = detector_it->second;

      if (iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
        ACTS_INFO ("Measurement is outlier on a fitter track, copy it over");
        tsosTypes[index] = Outlier;
        index++;
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
          index++;
          continue;
        }

        if (lastrot != nullptr && lastbutonerot != nullptr) {

          std::cout << lastrotindex << std::endl;
          std::cout << firstisshared << std::endl;

        }

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

        lastbutonerot = lastrot;
        lastrot       = &iMeasurement;
        lastrotindex  = index;   

        std::cout << lastrotindex << std::endl; 
        index++;
        continue;      
      }

      ACTS_VERBOSE ("Measurement is not shared, Reject it");
      tsosTypes[index] = RejectedHit;
      index++;
    }

    std::vector<std::size_t> newMeasurementsPerTrack;
    std::size_t measurment = 0;

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
      else {
        std::cout << "Try to recover shared hit" << std::endl;
        if (tracksPerMeasurement[measurment].size() > m_maxSharedTracksPerMeasurement){
          std::cout << "Too many tracks sharing this hit, drop it" << std::endl;
        }
        else {
          std::cout << "Keep shared hit" << std::endl;
          newMeasurementsPerTrack.push_back(measurment);
        }
      }
    }


    std::size_t sharedMeasurementsPerTrack = 0;
    for (auto iMeasurement : newMeasurementsPerTrack) {
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
        ++sharedMeasurementsPerTrack;
      }
    }
  
  

    if (sharedMeasurementsPerTrack > 0 && trackScore[iTrack] < m_minScoreSharedTracks) {
      TrkCouldBeAccepted = false;
      std::cout << "Track " << iTrack << " could not be accepted" << std::endl;
      iTrack++;
      continue;
    }
    
    for ( std::size_t i = 0; i< volumeList.size(); ++i){
      auto detector_it = m_volumeMap.find(volumeList[i]);
      if(detector_it == m_volumeMap.end()){
        continue;
      }
      auto detector = detector_it->second;
      if (counterMap[detector.detectorId].nSharedHits > detector.maxSharedHits){
        TrkCouldBeAccepted = false;
      }
    }

    if (TrkCouldBeAccepted){
      cleanTracks.push_back(iTrack);
      // ACTS_INFO("Track " << iTrack << " is clean");
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Track " << iTrack << " has " << sharedMeasurementsPerTrack << " shared measurements" << std::endl;
      std::cout << "Track " << iTrack << " has " << measurementsPerTrack[iTrack].size() << " measurements" << std::endl;
      std::cout << "Track " << iTrack << " has " << trackScore[iTrack] << " score" << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    iTrack++;
    continue;
    }
    iTrack++;
  }
  return cleanTracks;
}


