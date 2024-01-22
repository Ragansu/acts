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
      std::string name;
      std::size_t hits_score = 0;
      std::size_t holes_score = 0;
      std::size_t outliers_score = 0;
      std::size_t other_score = 0;
      std::size_t score = 0;
      unsigned int layerMin = 0;
      unsigned int layerMax = 0;

      DectectorConfig(std::string name, std::size_t hits_score, std::size_t holes_score, std::size_t outliers_score, std::size_t other_score, unsigned int layerMin, unsigned int layerMax) {
        this->name = name;
        this->hits_score = hits_score;
        this->holes_score = holes_score;
        this->outliers_score = outliers_score;
        this->other_score = other_score;
        this->layerMin = layerMin;
        this->layerMax = layerMax;
        this->score = 0;
      }

      void setVolumeIds(std::vector<unsigned int> volumeIds) {
        m_volumeIds = volumeIds;
      }
      void setLayerIds(std::vector<unsigned int> layerIds) {
        m_layerIds = layerIds;
      }
    protected:
      std::vector<unsigned int> getVolumeIds() {
        return m_volumeIds;
      std::vector<unsigned int> getLayerIds() {
        return m_layerIds;
      }
      
    private:
      std::vector<unsigned int> m_volumeIds = {};
      std::vector<unsigned int> m_layerIds = {};

  };

  std::vector<DectectorConfig> m_detectorConfigs = {
    DectectorConfig("Pixel", 20, -10, -5, 0, 0, 2),
    DectectorConfig("SCT", 10, -5, -2, 0, 0, 2),
    DectectorConfig("TRT", 0, 0, 0, 0, 0, 2),
    DectectorConfig("MDT", 20, 0, 0, 0, 0, 2),
    DectectorConfig("TGC", 20, 0, 0, 0, 0, 2),
    DectectorConfig("CSC", 20, 0, 0, 0, 0, 2),
    DectectorConfig("RPC", 20, 0, 0, 0, 0, 2)
  };

  class Detector : public DectectorConfig {
    std::size_t nMeasurements = 0;
    std::size_t nOutliers = 0;
    std::size_t nHoles = 0;

    std::vector<unsigned int> measurementVolume = {};
    std::vector<unsigned int> measurementLayer = {};
    std::vector<unsigned int> outlierVolume = {};
    std::vector<unsigned int> outlierLayer = {};
    std::vector<unsigned int> holeVolume = {};
    std::vector<unsigned int> holeLayer = {};

  };
  
  std::vector<Detector> m_detectors;
  };
  std::vector<Detector> constructDetectors(const ConstTrackContainer& tracks) const;

  int m_nDetectors = m_detectors.size();

  
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

  


  TypeScore m_typeScore;
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

  void fillTrajectoryState(const ConstTrackContainer& tracks, std::size_t trackID) const;
  int simpleScore(const ConstTrackContainer& tracks) const;
  void getCleanedOutTracks(const ConstTrackContainer& tracks) const;
  std::vector<std::size_t> solveAmbiguity(const ConstTrackContainer& tracks, int trackScore) const;
};

}  // namespace ActsExamples
