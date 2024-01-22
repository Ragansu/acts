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

  class DectectorConfig {
    public:
      DectectorConfig( std::size_t hits_score, std::size_t holes_score, std::size_t outliers_score, std::size_t other_score) {
        this->hits_score = hits_score;
        this->holes_score = holes_score;
        this->outliers_score = outliers_score;
        this->other_score = other_score;
      }
      
      std::size_t getHitsScore() const { return hits_score; }
      std::size_t getHolesScore() const { return holes_score; }
      std::size_t getOutliersScore() const { return outliers_score; }
      std::size_t getOtherScore() const { return other_score; }

    private:
      std::size_t hits_score;
      std::size_t holes_score;
      std::size_t outliers_score;
      std::size_t other_score;
  };

  std::map<unsigned int,DectectorConfig> Volumemap {
    {0,DectectorConfig( 20, -10, -2, 0)}, //pixel
    {1,DectectorConfig( 10, -5, -2, 0)},  //sct

  };

// muons TODO etahits and phihits

  std::size_t m_minScore = 0;
  
  // struct TypeScore {
  //   int value;
  //   std::string name;
  // };

  // std::vector<TypeScore> m_typeScores = {
  //   {20,"numberOfPixelHits"},
  //   {-10,"numberOfPixelHoles"},
  //   {10,"numberOfInnermostPixelLayerHits"},
  //   {-5,"numberOfGangedPixels"},
  //   {10,"numberOfSCTHits"},
  //   {-5,"numberOfSCTHoles"},
  //   {-2,"numberOfOutliersOnTrack"},
  //   {20,"numberOfMdtHits"},
  //   {20,"numberOfTgcPhiHits"},
  //   {10,"numberOfTgcEtaHits"},
  //   {20,"numberOfCscPhiHits"},
  //   {20,"numberOfCscEtaHits"},
  //   {20,"numberOfRpcPhiHits"},
  //   {10,"numberOfRpcEtaHits"}
  // };


 protected:
  /// Associated measurements ID to Tracks ID
  ///
  /// @param tracks is the input track container
  /// @param nMeasurementsMin minimum number of measurement per track
  /// @return an ordered list containing pairs of track ID and associated measurement ID
  // std::multimap<int, std::pair<std::size_t, std::vector<std::size_t>>>
  // mapTrackHits(const ConstTrackContainer& tracks, int nMeasurementsMin) const;

  /// Prepare the output track container to be written
  ///
  /// @param tracks is the input track container
  /// @param goodTracks is list of the IDs of all the tracks we want to keep
  ConstTrackContainer prepareOutputTrack(
      const ConstTrackContainer& tracks,
      std::vector<std::size_t>& goodTracks) const;

  std::vector<int> simpleScore(const ConstTrackContainer& tracks) const;
  std::vector<std::size_t> getCleanedOutTracks(const ConstTrackContainer& tracks) const;
  std::vector<std::size_t> solveAmbiguity(const ConstTrackContainer& tracks, std::vector<int> trackScore) const;
};

}  // namespace ActsExamples
