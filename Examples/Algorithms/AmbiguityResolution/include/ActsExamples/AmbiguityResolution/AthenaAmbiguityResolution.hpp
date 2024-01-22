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
      std::string name; // name of the detector
      std::size_t hits_score = 0; // score for hits
      std::size_t holes_score = 0; // score for holes
      std::size_t outliers_score = 0; // score for outliers
      std::size_t eta_score = 0; // score based on eta
      std::size_t phi_score = 0; // score based on phi
      std::size_t other_score = 0; // score for other measurements

      unsigned int layerMin = 0;
      unsigned int layerMax = 0;

      DectectorConfig(std::string name, std::size_t hits_score, std::size_t holes_score, std::size_t outliers_score, std::size_t eta_score, std::size_t phi_score, std::size_t other_score) 
      : name(name), hits_score(hits_score), holes_score(holes_score), outliers_score(outliers_score), eta_score(eta_score), phi_score(phi_score), other_score(other_score) { }
      void addAllTracks(std::size_t trackID) {
        AllTracks.push_back(trackID);
      }

      void addGoodTracks(std::size_t trackID) {
        GoodTracks.push_back(trackID);
      }

      // void setVolumeIds(std::vector<unsigned int> volumeIds) {
      //   m_volumeIds = volumeIds;
      // }
      // void setLayerIds(std::vector<unsigned int> layerIds) {
      //   m_layerIds = layerIds;
      }
    protected:
      // std::vector<unsigned int> getVolumeIds() {
      //   return m_volumeIds;
      // std::vector<unsigned int> getLayerIds() {
      //   return m_layerIds;
      // }
      
    private:
      // std::vector<unsigned int> m_volumeIds = {};
      // std::vector<unsigned int> m_layerIds = {};

      std::vector<std::size_t> AllTracks = {};
      std::vector<std::size_t> GoodTracks = {};

  };

  std::vector<DectectorConfig> m_detectorConfigs = {
    {"Pixel", 20, -10, -2, 0, 0, 0},
    {"SCT", 10, -5, -2, 0, 0, 0},
    {"MDT", 20,0,0,0,0,0},
    {"TGC", 0,0,0,20,10,0},
    {"CSC", 0,0,0,20,10,0},
    {"RPC", 0,0,0,20,10,0}
  };

  class Detector : public DectectorConfig {
  public:
    std::size_t tipIndex = 0;
    std::size_t nMeasurements = 0;
    std::size_t nOutliers = 0;
    std::size_t nHoles = 0;
    std::size_t score = 0;

    std::vector<unsigned int> measurementVolume = {};
    std::vector<unsigned int> measurementLayer = {};
    std::vector<unsigned int> outlierVolume = {};
    std::vector<unsigned int> outlierLayer = {};
    std::vector<unsigned int> holeVolume = {};
    std::vector<unsigned int> holeLayer = {};

    Detector(const DectectorConfig other) : DectectorConfig(other) { }

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

  std::vector<int> simpleScore(const ConstTrackContainer& tracks) const;
  void getCleanedOutTracks(const ConstTrackContainer& tracks) const;
  std::vector<std::size_t> solveAmbiguity(const ConstTrackContainer& tracks, int trackScore) const;
};

}  // namespace ActsExamples
