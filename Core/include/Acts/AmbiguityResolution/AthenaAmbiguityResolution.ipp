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

  std::cout << "Starting to compute initial state" << std::endl;

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

    int score = 100;
    auto counterMap = m_counterMap;
    
    if (track.chi2() > 0 && track.nDoF() > 0) {
      score+= 10*(1.0-std::log10(track.chi2()/track.nDoF())); // place holder
    }

    // TODO: add scored based on chi2 and ndof

    // detector score is determined by the number of hits/hole/outliers * hit/hole/outlier score
    // here so instead of calculating nHits/nHoles/nOutliers per volume, 
    // we just loop over all volumes and add the score.

    for (auto ts : track.trackStatesReversed()) {
      auto iVolume = ts.referenceSurface().geometryId().volume();
      auto iTypeFlags = ts.typeFlags();  

      auto detector_it = m_cfg.volumeMap.find(iVolume);
      if(detector_it != m_cfg.volumeMap.end()){
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
      else {
        ACTS_INFO("Detector not found");
        ACTS_INFO("Volume: " << iVolume);
      }
    }
    counterMaps.push_back(counterMap);
        if (Acts::VectorHelpers::phi(track.momentum()) > m_cfg.phiMax || 
        Acts::VectorHelpers::phi(track.momentum()) < m_cfg.phiMin) {
      score = 0;
    }

    if (Acts::VectorHelpers::eta(track.momentum()) > m_cfg.etaMax || 
        Acts::VectorHelpers::eta(track.momentum()) < m_cfg.etaMin) {
      score = 0;
    }
    iTrack++;
    trackScore.push_back(score);
    ACTS_INFO("Track score: " << score);

  } // end of loop over tracks
    
  return trackScore;
}


template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int>
Acts::AthenaAmbiguityResolution::solveAmbiguity(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack) const {

  std::cout << "Solving ambiguity" << std::endl;
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