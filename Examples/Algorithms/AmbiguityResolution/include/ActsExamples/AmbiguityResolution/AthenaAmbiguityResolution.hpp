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
  struct Detector {
    std::string name;
    std::size_t nMeasurements = 0;
    std::size_t nOutliers = 0;
    std::size_t nHoles = 0;

    std::vector<unsigned int> measurementVolume = {};
    std::vector<unsigned int> measurementLayer = {};
    std::vector<unsigned int> outlierVolume = {};
    std::vector<unsigned int> outlierLayer = {};
    std::vector<unsigned int> holeVolume = {};
    std::vector<unsigned int> holeLayer = {};

    void findnHoles(const ConstTrackContainer& tracks, std::size_t trackID) const;
    void findnOutliers(const ConstTrackContainer& tracks, std::size_t trackID) const;
    void findnMeasurements(const ConstTrackContainer& tracks, std::size_t trackID) const;
  };
  
  std::vector<Detector> m_detectors = {
    {"Pixel",0,0,0,{},{}},
    {"SCT",0,0,0,{},{}},
    {"MDT",0,0,0,{},{}},
    {"TGC",0,0,0,{},{}},
    {"CSC",0,0,0,{},{}},
    {"RPC",0,0,0,{},{}}

  };
  
  struct TypeScore {
    int value;
    std::string name;
  };

  std::vector<TypeScore> m_typeScores = {
    {20,"numberOfPixelHits"},
    {-10,"numberOfPixelHoles"},
    {10,"numberOfInnermostPixelLayerHits"},
    {-5,"numberOfGangedPixels"},
    {10,"numberOfSCTHits"},
    {-5,"numberOfSCTHoles"},
    {-2,"numberOfOutliersOnTrack"},
    {20,"numberOfMdtHits"},
    {20,"numberOfTgcPhiHits"},
    {10,"numberOfTgcEtaHits"},
    {20,"numberOfCscPhiHits"},
    {20,"numberOfCscEtaHits"},
    {20,"numberOfRpcPhiHits"},
    {10,"numberOfRpcEtaHits"}
  };

  


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
