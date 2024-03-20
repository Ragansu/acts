// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"


#include <unordered_map>

namespace Acts {

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
const TrackContainer<track_container_t, traj_t, holder_t>
Acts::AthenaAmbiguityResolution::prepareOutputTrack(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::size_t>& goodTracks) const {
  auto trackStateContainer =
      tracks.trackStateContainerHolder();
  auto trackContainer = std::make_shared<VectorTrackContainer>();
  trackContainer->reserve(goodTracks.size());
  // Temporary empty track state container: we don't change the original one,
  // but we need one for filtering
  auto tempTrackStateContainer =
      std::make_shared<VectorMultiTrajectory>();

  TrackContainer solvedTracks{trackContainer, tempTrackStateContainer};
  solvedTracks.ensureDynamicColumns(tracks);

  for (auto&& iTrack : goodTracks) {
    auto destProxy = solvedTracks.getTrack(solvedTracks.addTrack());
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFrom(srcProxy, false);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  const TrackContainer<track_container_t, traj_t, holder_t> outputTracks{
      std::make_shared<VectorTrackContainer>(
          std::move(*trackContainer)),
      trackStateContainer};
  return outputTracks;
}


template <typename track_container_t, typename traj_t,
          template <typename> class holder_t, typename source_link_hash_t,
          typename source_link_equality_t>
std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> AthenaAmbiguityResolution::computeInitialState(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    source_link_hash_t&& sourceLinkHash,
    source_link_equality_t&& sourceLinkEquality) const {
  auto measurementIndexMap =
      std::unordered_map<SourceLink, std::size_t, source_link_hash_t,
                         source_link_equality_t>(0, sourceLinkHash,
                                                 sourceLinkEquality);

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack;

  ACTS_INFO ( "Starting to compute initial state" << std::endl);

  int numberOfTracks = 0;
  for (const auto& track : tracks) {

    std::vector<std::tuple<std::size_t, std::size_t, bool >> measurements_tuples;

    for (auto ts : track.trackStatesReversed()) {
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();

        const auto& geoID = ts.referenceSurface().geometryId();   

        // assign a new measurement index if the source link was not seen yet
        auto emplace = measurementIndexMap.try_emplace(
            sourceLink, measurementIndexMap.size());

        bool isoutliner = ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag);


        
        measurements_tuples.push_back(std::make_tuple(emplace.first->second, geoID.volume(), isoutliner));
      }
    }

    measurementsPerTrack.push_back(std::move(measurements_tuples));

    ++numberOfTracks;
  }

  return measurementsPerTrack;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int> Acts::AthenaAmbiguityResolution::simpleScore(
 const TrackContainer<track_container_t, traj_t, holder_t>& tracks, 
 std::vector<std::map<std::size_t, Counter>>& counterMaps) const {

  std::vector<int> trackScore;
  int iTrack = 0;  

  // Loop over all the trajectories in the events
  for (const auto& track : tracks){

    auto counterMap = m_counterMap;
    
    // detector score is determined by the number of hits/hole/outliers * hit/hole/outlier score
    // here so instead of calculating nHits/nHoles/nOutliers per volume, 
    // we just loop over all volumes and add the score.

    std::map<std::size_t,VolumeConfig> detectorMap;
    std::vector<std::size_t> detectorList;

    for (auto ts : track.trackStatesReversed()) {
      auto iVolume = ts.referenceSurface().geometryId().volume();
      auto iTypeFlags = ts.typeFlags();  

      auto detector_it = m_cfg.volumeMap.find(iVolume);
      if(detector_it != m_cfg.volumeMap.end()){
        auto detector = detector_it->second;

        auto volume_it = detectorMap.find(detector.detectorId);
        if(volume_it == detectorMap.end()){
          detectorMap[detector.detectorId] = detector;
          detectorList.push_back(detector.detectorId);
        }

        if (iTypeFlags.test(Acts::TrackStateFlag::MeasurementFlag)){
          if (iTypeFlags.test(Acts::TrackStateFlag::SharedHitFlag)){
            counterMap[detector.detectorId].nSharedHits++;
          }
          counterMap[detector.detectorId].nHits++;
        }
        else if (iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
          counterMap[detector.detectorId].nOutliers++;
        }
        else if (iTypeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
          counterMap[detector.detectorId].nHoles++;
        }
      }
      else {
        ACTS_WARNING("Detector not found at Volume: " << iVolume);
      }
    }
    counterMaps.push_back(counterMap);

    ACTS_VERBOSE ( "Number of detectors: " << detectorList.size() << std::endl);

    bool TrkCouldBeAccepted = true;
    int score = 0;

    if (!TrkCouldBeAccepted){
      ACTS_VERBOSE ( "Track " << iTrack << " could not be accepted" << std::endl);
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_VERBOSE("Track: "<< iTrack <<" score: " << score << " could not be accepted");
      continue;    
    }

    if (Acts::VectorHelpers::perp(track.momentum()) > m_cfg.pTMax || 
        Acts::VectorHelpers::perp(track.momentum()) < m_cfg.pTMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_VERBOSE("Track: "<< iTrack <<" score: " << score << " pT: " << Acts::VectorHelpers::perp(track.momentum()));
      continue;
    } 

    if (Acts::VectorHelpers::phi(track.momentum()) > m_cfg.phiMax || 
        Acts::VectorHelpers::phi(track.momentum()) < m_cfg.phiMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_VERBOSE("Track: "<< iTrack <<" score: " << score << " phi: " << Acts::VectorHelpers::phi(track.momentum()));
      continue;
    }

    if (Acts::VectorHelpers::eta(track.momentum()) > m_cfg.etaMax || 
        Acts::VectorHelpers::eta(track.momentum()) < m_cfg.etaMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_VERBOSE("Track: "<< iTrack <<" score: " << score << " eta: " << Acts::VectorHelpers::eta(track.momentum()));
      continue;   
    }

    // real scoring starts here
    if (m_useAmbigFcn) {
      score = 0;
    }
    else {

      score = 100;

      for(long unsigned int i = 0; i< detectorList.size(); i++){

        auto volume_it = detectorMap.find(detectorList[i]);
        auto detector = volume_it->second;
        
        ACTS_DEBUG ("---> Found summary information");
        ACTS_DEBUG("---> Detector ID: " << detector.detectorId);
        ACTS_DEBUG ("---> Number of hits: " << counterMap[detector.detectorId].nHits);
        ACTS_DEBUG ("---> hitsScoreWeight: " << detector.hitsScoreWeight);
        ACTS_DEBUG ("---> Number of holes: " << counterMap[detector.detectorId].nHoles);
        ACTS_DEBUG ("---> holesScoreWeight: " << detector.holesScoreWeight); 
        ACTS_DEBUG ("---> Number of outliers: " << counterMap[detector.detectorId].nOutliers);
        ACTS_DEBUG ("---> outliersScoreWeight: " << detector.outliersScoreWeight);
        ACTS_DEBUG ("---> otherScoreWeight: " << detector.otherScoreWeight);

        ACTS_DEBUG ("---> Min hits: " << detector.minHits);
        ACTS_DEBUG ("---> Max holes: " << detector.maxHoles);
        ACTS_DEBUG ("---> Max outliers: " << detector.maxOutliers);
        ACTS_DEBUG ("---> Max unused: " << detector.maxUnused);

        if (counterMap[detector.detectorId].nHits < detector.minHits){
          TrkCouldBeAccepted = false;
        }

        else if (counterMap[detector.detectorId].nHoles > detector.maxHoles){
          TrkCouldBeAccepted = false;
        }

        else if (counterMap[detector.detectorId].nOutliers > detector.maxOutliers){
          TrkCouldBeAccepted = false;
        }
        else {

          score += counterMap[detector.detectorId].nHits * detector.hitsScoreWeight;
          score += counterMap[detector.detectorId].nHoles * detector.holesScoreWeight;
          score += counterMap[detector.detectorId].nOutliers * detector.outliersScoreWeight;
          score += counterMap[detector.detectorId].nSharedHits * detector.otherScoreWeight;
        }
      }

      if (track.chi2() > 0 && track.nDoF() > 0) {
        
        double p = 1. / log10 (10. + 10. * track.chi2() / track.nDoF());
        if (p > 0) {
          score+= p;
        }
        else score -= 50;
      } 
    }   

    iTrack++;
    trackScore.push_back(score);
    ACTS_VERBOSE("Track: "<< iTrack <<" score: " << score);

  } // end of loop over tracks

  if (m_useAmbigFcn) {
    ACTS_VERBOSE ( "Using ambiguity function" << std::endl);
    trackScore = ambigScore(tracks, trackScore);
  }
  else {
    ACTS_VERBOSE ( "Not using ambiguity function" << std::endl);
  }
    
  return trackScore;
}


template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int> Acts::AthenaAmbiguityResolution::ambigScore(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<int> trackScore) const {

  std::vector<int> trackScoreAmbig;
  int iTrack = 0;
  for (auto track : tracks) {
    double score = trackScore[iTrack];
    if (score >= 0) {
      trackScoreAmbig.push_back(score);
      iTrack++;
      continue;
    }

    double pT = Acts::VectorHelpers::perp(track.momentum());

    double prob = log10(pT);

    if (track.chi2() > 0 && track.nDoF() > 0) {
      double p = 1. / log10 (10. + 10. * track.chi2() / track.nDoF());
      prob *= p;
    }
   
    trackScoreAmbig.push_back(prob);
  }  // work in progress
  return trackScoreAmbig;
}


template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int>
Acts::AthenaAmbiguityResolution::solveAmbiguity(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack) const {

  ACTS_INFO ( "Solving ambiguity" << std::endl);  
  ACTS_INFO ( "Number of tracks: " << tracks.size() << std::endl);
  ACTS_INFO ( "Config file location: " << m_cfg.volumeFile << std::endl);
  std::vector<std::map<std::size_t, Counter>> counterMaps;
  std::vector<int> trackScore = simpleScore(tracks, counterMaps);

  std::vector<std::size_t> cleanTracks = getCleanedOutTracks(trackScore, counterMaps, measurementsPerTrack);

  ACTS_INFO("Number of tracks: " << tracks.size());
  ACTS_INFO("Number of clean tracks: " << cleanTracks.size());
  ACTS_INFO("Min score: " << m_cfg.minScore);

  std::vector<int> goodTracks;
  std::size_t iTrack = 0;
  for (auto track : tracks) {
    if (trackScore[iTrack] > m_cfg.minScore) {
      goodTracks.push_back(track.index());
    }
    iTrack++;
  }

  ACTS_INFO("Number of good tracks: " << goodTracks.size());
  return goodTracks;
}

}  // namespace Acts