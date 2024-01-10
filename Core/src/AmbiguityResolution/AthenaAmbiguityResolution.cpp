// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"

namespace Acts {

namespace {

/// Removes a track from the state which has to be done for multiple properties
/// because of redundancy.
static void removeTrack(AthenaAmbiguityResolution::State& state,
                        std::size_t iTrack) {
  for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
    state.tracksPerMeasurement[iMeasurement].erase(iTrack);

    if (state.tracksPerMeasurement[iMeasurement].size() == 1) {
      auto jTrack = *state.tracksPerMeasurement[iMeasurement].begin();
      --state.sharedMeasurementsPerTrack[jTrack];
    }
  }

  state.selectedTracks.erase(iTrack);
}

}  // namespace

void AthenaAmbiguityResolution::resolve(State& state) const {
  /// Compares two tracks based on the number of shared measurements in order to
  /// decide if we already met the final state.
  auto sharedMeasurementsComperator = [&state](std::size_t a, std::size_t b) {
    return state.sharedMeasurementsPerTrack[a] <
           state.sharedMeasurementsPerTrack[b];
  };

  /// Compares two tracks in order to find the one which should be evicted.
  /// First we compare the relative amount of shared measurements. If that is
  /// indecisive we use the chi2.
  auto trackComperator = [&state](std::size_t a, std::size_t b) {
    /// Helper to calculate the relative amount of shared measurements.
    auto relativeSharedMeasurements = [&state](std::size_t i) {
      return 1.0 * state.sharedMeasurementsPerTrack[i] /
             state.measurementsPerTrack[i].size();
    };

    if (relativeSharedMeasurements(a) != relativeSharedMeasurements(b)) {
      return relativeSharedMeasurements(a) < relativeSharedMeasurements(b);
    }
    return state.trackChi2[a] < state.trackChi2[b];
  };

  for (std::size_t i = 0; i < m_cfg.maximumIterations; ++i) {
    // Lazy out if there is nothing to filter on.
    if (state.selectedTracks.empty()) {
      ACTS_VERBOSE("no tracks left - exit loop");
      break;
    }

    // Find the maximum amount of shared measurements per track to decide if we
    // are done or not.
    auto maximumSharedMeasurements = *std::max_element(
        state.selectedTracks.begin(), state.selectedTracks.end(),
        sharedMeasurementsComperator);
    ACTS_VERBOSE(
        "maximum shared measurements "
        << state.sharedMeasurementsPerTrack[maximumSharedMeasurements]);
    if (state.sharedMeasurementsPerTrack[maximumSharedMeasurements] <
        m_cfg.maximumSharedHits) {
      break;
    }

    // Find the "worst" track by comparing them to each other
    auto badTrack =
        *std::max_element(state.selectedTracks.begin(),
                          state.selectedTracks.end(), trackComperator);
    ACTS_VERBOSE("remove track "
                 << badTrack << " nMeas "
                 << state.measurementsPerTrack[badTrack].size() << " nShared "
                 << state.sharedMeasurementsPerTrack[badTrack] << " chi2 "
                 << state.trackChi2[badTrack]);
    removeTrack(state, badTrack);
  }
}

void AthenaAmbiguityResolution::trackScoringTool() const {
  constexpr int toZero{ 0 };
  constexpr size_t numberOfInDetCounters{ 14 };
  constexpr std::array<size_t, numberOfInDetCounters> inDetIndex{
    numberOfPixelHits,
    numberOfPixelHoles,
    numberOfInnermostPixelLayerHits,
    numberOfGangedPixels,
    numberOfSCTHits,
    numberOfSCTHoles,
    numberOfOutliersOnTrack
    numberOfMdtHits,          
    numberOfTgcPhiHits,
    numberOfTgcEtaHits,       
    numberOfCscPhiHits,
    numberOfCscEtaHits,       
    numberOfRpcPhiHits,       
    numberOfRpcEtaHits,

  };
  setTheseElements(inDetIndex, toZero);
  
  //set some test values
	m_summaryTypeScore[numberOfPixelHits]	=  20;
	m_summaryTypeScore[numberOfPixelHoles] = -10;  // a hole is bad
	m_summaryTypeScore[numberOfInnermostPixelLayerHits] =  10;  // addition for being b-layer
	m_summaryTypeScore[numberOfGangedPixels] =  -5;  // decrease for being ganged
	m_summaryTypeScore[numberOfSCTHits] =  10;  // half of a pixel, since only 1dim
	m_summaryTypeScore[numberOfSCTHoles] =  -5;  // a hole is bad !
	m_summaryTypeScore[numberOfOutliersOnTrack] =  -2;  // an outlier might happen

	// scoring for Muons is missing
	m_summaryTypeScore[numberOfMdtHits]	= 20;
	m_summaryTypeScore[numberOfTgcPhiHits]	= 20;
	m_summaryTypeScore[numberOfTgcEtaHits]	= 10;
	m_summaryTypeScore[numberOfCscPhiHits]	= 20;
	m_summaryTypeScore[numberOfCscEtaHits]	= 20;
	m_summaryTypeScore[numberOfRpcPhiHits]	= 20;
	m_summaryTypeScore[numberOfRpcEtaHits]	= 10;
}





void AthenaAmbiguityResolution::simpleScore(State& state) const {
  TrackScore score = 0;
  for(j=0;j<state.numberOfTracks;j++) {
    state.trackScore[j] = 100; // score of 100 per track

    if (track.fitQuality() && state.trackDOF[j] < 0) {
      ACTS_VERBOSE("numberDoF < 0, reject it");
      state.trackScore[j] = 0;
      break;
    }

    // --- prob(chi2,NDF), protect for chi2<0
    if (track.fitQuality()!=nullptr && state.trackChi2[j] > 0 && state.trackDOF[j] > 0) {
      state.trackScore[j]+= std::log10(1.0-Genfun::CumulativeChiSquare(state.trackDOF[j])(state.trackChi2[j]));
    }
    int numberOfTrackSummaryTypes = 14;
    // --- summary score analysis
    for (int i=0; i<numberOfTrackSummaryTypes; ++i) {
      int value = trackSummary.get(static_cast<Trk::SummaryType>(i));
      //value is -1 if undefined.
      if (value>0) {
        state.trackScore[j]+=m_summaryTypeScore[i]*value;
        ACTS_VERBOSE("\tType ["<<i<<"], value \t= "<<value<<"], score \t="<<score);
      }
    }
    score+= state.trackScore[j];
  } // end of loop over tracks
  state.score = score;
}



// namespace Acts
