// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// #include "Acts/EventData/Track.hpp"
// #include "Acts/Framework/IAlgorithm.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <string>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

namespace Acts {

/// Generic implementation of the machine learning ambiguity resolution
/// Contains method for data preparations
class AthenaAmbiguityResolution {
 public:
  /// Framework execute method of the algorithm

  struct VolumeConfig {

      int hitsScoreWeight;
      int holesScoreWeight;
      int outliersScoreWeight;
      int otherScoreWeight;

      int minHits;
      int maxHoles;
      int maxOutliers;
      int maxUnused;
      int maxSharedHits;

      bool sharedHitsFlag;
      
      std::size_t detectorId;
  };

  struct Config {
  
    std::map<unsigned int,VolumeConfig> volumeMap;

    int minScore = 0;

    int minScoreSharedTracks = 0;

    std::size_t maxSharedTracksPerMeasurement = 10;

    std::size_t maxShared = 5;

    double phiMin = -M_PI;
    double phiMax = M_PI;

    double etaMin = -5;
    double etaMax = 5;
  };


  struct Counter {
    int nHits;
    int nHoles;
    int nOutliers;
    int nUnused;

    int nSharedHits;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param name name of the algorithm
  /// @param lvl is the logging level
  AthenaAmbiguityResolution(const Config& cfg,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("AthenaAmbiguityResolution",
                                                 Logging::INFO))
      : m_cfg{cfg}, m_logger{std::move(logger)} {}

// muons TODO etahits and phihits
  
template <typename track_container_t, typename traj_t,
          template <typename> class holder_t, typename source_link_hash_t,
          typename source_link_equality_t>
std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> computeInitialState(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    source_link_hash_t&& sourceLinkHash,
    source_link_equality_t&& sourceLinkEquality) const;



  /// Prepare the output track container to be written
  ///
  /// @param tracks is the input track container
  /// @param goodTracks is list of the IDs of all the tracks we want to keep
template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
const TrackContainer<track_container_t, traj_t, holder_t> prepareOutputTrack(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::size_t>& goodTracks) const;

  /// Compute the score of each track
  ///
  /// @param tracks is the input track container
  /// @return a vector of scores for each track
template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int> simpleScore(
 const TrackContainer<track_container_t, traj_t, holder_t>& tracks, 
 std::vector<std::map<std::size_t, Counter>>& counterMaps) const;
  /// Remove tracks that are not good enough based on cuts
  ///
  /// @param trackScore is the score of each track
  /// @param counterMaps is the counter map for each track
  /// @param measurementsPerTrack is the list of measurements for each track
  /// @return a vector of IDs of the tracks we want to keep
  std::vector<std::size_t> getCleanedOutTracks(
    std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack) const;

  /// Remove tracks that are not good enough
  ///
  /// @param tracks is the input track container
  /// @param trackScore is the score of each track
  /// @return a vector of IDs of the tracks we want to keep
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t>
  std::vector<int> solveAmbiguity(
      const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
      std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack) const;

private:
  Config m_cfg;
  
  std::map<std::size_t, Counter> m_counterMap = {
    {0,{0,0,0,0,0}},
    {1,{0,0,0,0,0}},
    {2,{0,0,0,0,0}}
  };
  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.ipp"