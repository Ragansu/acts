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
  struct TrackSummary {
    std::vector< uint8_t >         numberOfContribPixelLayers       ;
    std::vector< uint8_t >         numberOfBLayerHits               ;
    std::vector< uint8_t >         numberOfBLayerOutliers           ;
    std::vector< uint8_t >         numberOfBLayerSharedHits         ;
    std::vector< uint8_t >         numberOfBLayerSplitHits          ;
    std::vector< uint8_t >         expectBLayerHit                  ; /// @todo FIXME!Should be bool.
    std::vector< uint8_t >         numberOfInnermostPixelLayerHits               ;
    std::vector< uint8_t >         numberOfInnermostPixelLayerOutliers           ;
    std::vector< uint8_t >         numberOfInnermostPixelLayerSharedHits         ;
    std::vector< uint8_t >         numberOfInnermostPixelLayerSplitHits          ;
    std::vector< uint8_t >         expectInnermostPixelLayerHit                  ; /// @todo FIXME!Should be bool.
    std::vector< uint8_t >         numberOfNextToInnermostPixelLayerHits               ;
    std::vector< uint8_t >         numberOfNextToInnermostPixelLayerOutliers           ;
    std::vector< uint8_t >         numberOfNextToInnermostPixelLayerSharedHits         ;
    std::vector< uint8_t >         numberOfNextToInnermostPixelLayerSplitHits          ;
    std::vector< uint8_t >         expectNextToInnermostPixelLayerHit                  ; /// @todo FIXME!Should be bool.
    std::vector< uint8_t >         numberOfPixelHits                ;
    std::vector< uint8_t >         numberOfPixelOutliers            ;
    std::vector< uint8_t >         numberOfPixelHoles               ;
    std::vector< uint8_t >         numberOfPixelSharedHits          ;
    std::vector< uint8_t >         numberOfPixelSplitHits           ;
    std::vector< uint8_t >         numberOfGangedPixels             ;
    std::vector< uint8_t >         numberOfGangedFlaggedFakes       ;
    std::vector< uint8_t >         numberOfPixelDeadSensors         ;
    std::vector< uint8_t >         numberOfPixelSpoiltHits          ;
    std::vector< uint8_t >         numberOfSCTHits                  ;
    std::vector< uint8_t >         numberOfSCTOutliers              ;
    std::vector< uint8_t >         numberOfSCTHoles                 ;
    std::vector< uint8_t >         numberOfSCTDoubleHoles           ;
    std::vector< uint8_t >         numberOfSCTSharedHits            ;
    std::vector< uint8_t >         numberOfSCTDeadSensors           ;
    std::vector< uint8_t >         numberOfSCTSpoiltHits            ;
    std::vector< uint8_t >         numberOfTRTHits                  ;
    std::vector< uint8_t >         numberOfTRTOutliers              ;
    std::vector< uint8_t >         numberOfTRTHoles                 ;
    std::vector< uint8_t >         numberOfTRTHighThresholdHits     ;
    std::vector< uint8_t >         numberOfTRTHighThresholdOutliers ;
    std::vector< uint8_t >         numberOfTRTDeadStraws            ;
    std::vector< uint8_t >         numberOfTRTTubeHits              ;
    std::vector< uint8_t >         numberOfTRTXenonHits             ;
    
    std::vector< uint8_t >         numberOfPrecisionLayers;
    std::vector< uint8_t >         numberOfPrecisionHoleLayers;
    std::vector< uint8_t >         numberOfPhiLayers;
    std::vector< uint8_t >         numberOfPhiHoleLayers;
    std::vector< uint8_t >         numberOfTriggerEtaLayers;
    std::vector< uint8_t >         numberOfTriggerEtaHoleLayers;
    
    std::vector< uint8_t >         numberOfOutliersOnTrack          ;
    std::vector< uint8_t >         standardDeviationOfChi2OS        ;
    std::vector< float >           eProbabilityComb;
    std::vector< float >           eProbabilityHT;
    std::vector< float >           eProbabilityToT;
    std::vector< float >           eProbabilityBrem;

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

  int simpleScore(const ConstTrackContainer& tracks) const;
  void getCleanedOutTracks(const ConstTrackContainer& tracks) const;
  std::vector<std::size_t> solveAmbiguity(const ConstTrackContainer& tracks, int trackScore) const;
};

}  // namespace ActsExamples
