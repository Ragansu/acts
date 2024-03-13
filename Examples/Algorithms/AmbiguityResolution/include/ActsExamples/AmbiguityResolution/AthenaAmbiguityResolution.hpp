// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <map>
#include <string>
#include <vector>

namespace ActsExamples {

/// Generic implementation of the machine learning ambiguity resolution
/// Contains method for data preparations
class AthenaAmbiguityResolution : public IAlgorithm {
 public:
  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param name name of the algorithm
  /// @param lvl is the logging level
  AthenaAmbiguityResolution(std::string name, Acts::Logging::Level lvl);
  ///

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

  struct Counter {
    int nHits;
    int nHoles;
    int nOutliers;
    int nUnused;

    int nSharedHits;
  };

// muons TODO etahits and phihits

  int m_minScore = 0;
  int m_minScoreSharedTracks = 200;
  std::size_t m_maxSharedTracksPerMeasurement = 3;
  std::size_t m_maxShared = 5;
  

 protected:
  template <typename source_link_hash_t,
          typename source_link_equality_t>
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, Acts::ConstTrackStateType>>> 
      computeInitialState(
        const ConstTrackContainer& tracks,
        source_link_hash_t&& sourceLinkHash,
        source_link_equality_t&& sourceLinkEquality) const;



  /// Prepare the output track container to be written
  ///
  /// @param tracks is the input track container
  /// @param goodTracks is list of the IDs of all the tracks we want to keep
  ConstTrackContainer prepareOutputTrack(
      const ConstTrackContainer& tracks,
      std::vector<std::size_t>& goodTracks) const;

  /// Compute the score of each track
  ///
  /// @param tracks is the input track container
  /// @return a vector of scores for each track
  std::vector<int> simpleScore(const ConstTrackContainer& tracks,   std::vector<std::map<std::size_t, Counter>>& counterMaps) const;

  /// Remove tracks that are not good enough based on cuts
  ///
  /// @param tracks is the input track container
  /// @return a vector of IDs of the tracks we want to keep
  std::vector<std::size_t> getCleanedOutTracks(const ActsExamples::ConstTrackContainer& tracks ,
  std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps) const;

  /// Remove tracks that are not good enough
  ///
  /// @param tracks is the input track container
  /// @param trackScore is the score of each track
  /// @return a vector of IDs of the tracks we want to keep
  std::vector<std::size_t> solveAmbiguity(const ConstTrackContainer& tracks, std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps) const;

private:
  // std::map<unsigned int,VolumeConfig> m_volumeMap {
  //   {16,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // pixel 1
  //   {17,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // pixel 2
  //   {18,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // pixel 3

  //   {23,{15, -8, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // short strip 1
  //   {24,{15, -8, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // short strip 2
  //   {25,{15, -8, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // short strip 3

  //   {28,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // long strip 1
  //   {29,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // long strip 2
  //   {30,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}} // long strip 3

  // };

  std::map<unsigned int,VolumeConfig> m_volumeMap {

    {8,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 1
    {9,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 2 (barrel)
    {10,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 3

    {13,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 1 
    {14,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 2
    {15,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 3
    {16,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 4 (barrel)
    {18,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 5
    {19,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 6
    {20,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 7

    {22,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 1
    {23,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 2 (barrel)
    {24,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 3

    {25,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // HGTD 1
    {2,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // HGTD 2



  };
  std::map<std::size_t, Counter> m_counterMap = {
    {0,{0,0,0,0,0}},
    {1,{0,0,0,0,0}},
    {2,{0,0,0,0,0}}
  };
};

}  // namespace ActsExamples

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolution.ipp"