// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <unordered_map>

namespace Acts {

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
const TrackContainer<track_container_t, traj_t, holder_t>
Acts::AthenaAmbiguityResolution::prepareOutputTrack(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::size_t>& goodTracks) const {
  auto trackStateContainer = tracks.trackStateContainerHolder();
  auto trackContainer = std::make_shared<VectorTrackContainer>();
  trackContainer->reserve(goodTracks.size());
  // Temporary empty track state container: we don't change the original one,
  // but we need one for filtering
  auto tempTrackStateContainer = std::make_shared<VectorMultiTrajectory>();

  TrackContainer solvedTracks{trackContainer, tempTrackStateContainer};
  solvedTracks.ensureDynamicColumns(tracks);

  for (auto&& iTrack : goodTracks) {
    auto destProxy = solvedTracks.getTrack(solvedTracks.addTrack());
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFrom(srcProxy, false);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  const TrackContainer<track_container_t, traj_t, holder_t> outputTracks{
      std::make_shared<VectorTrackContainer>(std::move(*trackContainer)),
      trackStateContainer};
  return outputTracks;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t, typename source_link_hash_t,
          typename source_link_equality_t>
std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
AthenaAmbiguityResolution::computeInitialState(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    source_link_hash_t&& sourceLinkHash,
    source_link_equality_t&& sourceLinkEquality) const {
  auto measurementIndexMap =
      std::unordered_map<SourceLink, std::size_t, source_link_hash_t,
                         source_link_equality_t>(0, sourceLinkHash,
                                                 sourceLinkEquality);

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
      measurementsPerTrack;

  ACTS_INFO("Starting to compute initial state");

  for (const auto& track : tracks) {
    std::vector<std::tuple<std::size_t, std::size_t, bool>> measurements_tuples;

    for (auto ts : track.trackStatesReversed()) {
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();

        const auto& geoID = ts.referenceSurface().geometryId();

        // assign a new measurement index if the source link was not seen yet
        auto emplace = measurementIndexMap.try_emplace(
            sourceLink, measurementIndexMap.size());

        bool isoutliner =
            ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag);

        measurements_tuples.push_back(
            std::make_tuple(emplace.first->second, geoID.volume(), isoutliner));
      }
    }

    measurementsPerTrack.push_back(std::move(measurements_tuples));
  }

  return measurementsPerTrack;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<double> Acts::AthenaAmbiguityResolution::simpleScore(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::map<std::size_t, Counter>>& counterMaps,
    Optional_cuts<track_container_t, traj_t, holder_t> optionalCuts) const {
  std::vector<double> trackScore;

  int iTrack = 0;

  ACTS_INFO("Number of detectors: " << m_cfg.detectorMap.size());

  ACTS_INFO("Starting to score tracks");

  // Loop over all the trajectories in the events
  for (const auto& track : tracks) {
    auto counterMap = std::map<std::size_t, Counter>();

    // detector score is determined by the number of hits/hole/outliers *
    // hit/hole/outlier score here so instead of calculating
    // nHits/nHoles/nOutliers per volume, we just loop over all volumes and add
    // the score.
    bool doubleFlag = false;

    for (const auto& ts : track.trackStatesReversed()) {
      auto iVolume = ts.referenceSurface().geometryId().volume();
      auto iTypeFlags = ts.typeFlags();

      auto volume_it = m_cfg.volumeMap.find(iVolume);
      if (volume_it != m_cfg.volumeMap.end()) {
        auto detectorId = volume_it->second;
        if (!iTypeFlags.test(Acts::TrackStateFlag::HoleFlag))
          doubleFlag = false;

        if (iTypeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
          if (doubleFlag) {
            counterMap[detectorId].nDoubleHoles++;
            doubleFlag = false;
          } else {
            doubleFlag = true;
          };
          counterMap[detectorId].nHoles++;
        } else if (iTypeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
          if (iTypeFlags.test(Acts::TrackStateFlag::SharedHitFlag)) {
            counterMap[detectorId].nSharedHits++;
          }
          counterMap[detectorId].nHits++;
        } else if (iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
          counterMap[detectorId].nOutliers++;
        }
      } else {
        ACTS_DEBUG("Detector not found at Volume: " << iVolume);
      }
    }
    counterMaps.push_back(counterMap);

    double score = 1;

    if (Acts::VectorHelpers::perp(track.momentum()) > m_cfg.pTMax ||
        Acts::VectorHelpers::perp(track.momentum()) < m_cfg.pTMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score: " << score << " pT: "
                           << Acts::VectorHelpers::perp(track.momentum()));
      continue;
    }

    if (Acts::VectorHelpers::phi(track.momentum()) > m_cfg.phiMax ||
        Acts::VectorHelpers::phi(track.momentum()) < m_cfg.phiMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score: " << score << " phi: "
                           << Acts::VectorHelpers::phi(track.momentum()));
      continue;
    }

    if (Acts::VectorHelpers::eta(track.momentum()) > m_cfg.etaMax ||
        Acts::VectorHelpers::eta(track.momentum()) < m_cfg.etaMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score: " << score << " eta: "
                           << Acts::VectorHelpers::eta(track.momentum()));
      continue;
    }

    for (const auto& ambicut : optionalCuts.cuts) {
      if (!ambicut(track)) {
        score = 0;
        break;
      }
    }

    if (score == 0) {
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score : " << score);
      continue;
    }

    for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
         detectorId++) {
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;

      ACTS_DEBUG("---> Found summary information");
      ACTS_DEBUG("---> Detector ID: " << detectorId);
      ACTS_DEBUG("---> Number of hits: " << counterMap[detectorId].nHits);
      ACTS_DEBUG("---> Number of holes: " << counterMap[detectorId].nHoles);
      ACTS_DEBUG("---> Number of double holes: "
                 << counterMap[detectorId].nDoubleHoles);
      ACTS_DEBUG(
          "---> Number of outliers: " << counterMap[detectorId].nOutliers);

      if ((counterMap[detectorId].nHits < detector.minHits) ||
          (counterMap[detectorId].nHits > detector.maxHits) ||
          (counterMap[detectorId].nHoles > detector.maxHoles) ||
          (counterMap[detectorId].nDoubleHoles > detector.maxDoubleHoles) ||
          (counterMap[detectorId].nOutliers > detector.maxOutliers)) {
        score = 0;
        ACTS_DEBUG("Track: " << iTrack << " score from cuts: " << score);
        break;
      }
    }

    if (score <= 0) {
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score : " << score);
      continue;
    }

    // real scoring starts here

    if (m_cfg.useAmbigFcn) {
      score = 1;
    } else {
      score = 100;
      for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
           detectorId++) {
        auto detector_it = m_cfg.detectorMap.find(detectorId);
        auto detector = detector_it->second;

      }

      for (const auto& ambiweights : optionalCuts.weights) {
        ambiweights(track, score);
      }

      if (track.chi2() > 0 && track.nDoF() > 0) {
        double p = 1. / log10(10. + 10. * track.chi2() / track.nDoF());
        if (p > 0) {
          score += p;
        } else
          score -= 50;
      }
    }

    iTrack++;
    trackScore.push_back(score);
    ACTS_VERBOSE("Track: " << iTrack << " score: " << score);

  }  // end of loop over tracks

  if (!m_cfg.useAmbigFcn) {
    ACTS_INFO("Not using ambiguity function");
    return trackScore;
  }

  ACTS_INFO("Using ambiguity function");

  std::vector<double> trackScoreAmbig;
  iTrack = 0;
  for (const auto& track : tracks) {
    double score = trackScore[iTrack];
    if (score == 0) {
      trackScoreAmbig.push_back(0.0f);
      iTrack++;
      continue;
    }
    auto counterMap = counterMaps[iTrack];
    double pT = Acts::VectorHelpers::perp(track.momentum());

    double prob = log10(pT * 1000) - 1.;
    // pT in GeV, hence 100 MeV is minimum and gets score = 1
    ACTS_DEBUG("Modifier for pT = " << pT << " GeV is : " << prob
                                    << "  New score now: " << prob);

    for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
         detectorId++) {
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;

      std::size_t iHits = counterMap[detectorId].nHits;
      if (detector.factorHits.size() < iHits) {
        ACTS_WARNING("Detector " << detectorId
                                 << " has not enough factors for hits");
        continue;
      }
      if (iHits < detector.minHits) {
        prob /= (detector.minHits - iHits + 1);  // missing hits are bad !
        iHits = detector.minHits;
      }
      prob *= detector.factorHits[iHits];
      ACTS_DEBUG("Modifier for " << iHits
                                 << " hits: " << detector.factorHits[iHits]
                                 << "  New score now: " << prob);

      std::size_t iHoles = counterMap[detectorId].nHoles;
      if (detector.factorHoles.size() < iHoles) {
        ACTS_WARNING("Detector " << detectorId
                                 << " has not enough factors for holes");
        continue;
      }
      if (iHoles > detector.maxHoles) {
        prob /= (iHoles - detector.maxHoles + 1);  // holes are bad !
        iHoles = detector.maxHoles;
      }
      prob *= detector.factorHoles[iHoles];
      ACTS_DEBUG("Modifier for " << iHoles
                                 << " holes: " << detector.factorHoles[iHoles]
                                 << "  New score now: " << prob);
    }

    for (const auto& ambiscore : optionalCuts.ambiscores) {
      ambiscore(track, prob);
    }

    if (track.chi2() > 0 && track.nDoF() > 0) {
      double chi2 = track.chi2();
      int indf = track.nDoF();
      double fac = 1. / log10(10. + 10. * chi2 / indf);
      prob *= fac;
      ACTS_DEBUG("Modifier for chi2 = " << chi2 << " and NDF = " << indf
                                        << " is : " << fac
                                        << "  New score now: " << prob)
    }
    trackScoreAmbig.push_back(prob);
  }  // work in progress
  return trackScoreAmbig;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int> Acts::AthenaAmbiguityResolution::solveAmbiguity(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
        measurementsPerTrack,
    Optional_cuts<track_container_t, traj_t, holder_t> optionalCuts) const {
  ACTS_INFO("Solving ambiguity");
  ACTS_INFO("Number of tracks: " << tracks.size());
  ACTS_INFO("Config file location: " << m_cfg.configFile);
  std::vector<std::map<std::size_t, Counter>> counterMaps;
  std::vector<double> trackScore =
      simpleScore(tracks, counterMaps, optionalCuts);

  for (const auto& track : tracks) {
    ACTS_INFO(
        "Track: " << track.index()
                  << " pT: " << Acts::VectorHelpers::perp(track.momentum())
                  << " eta: " << Acts::VectorHelpers::eta(track.momentum())
                  << " phi: " << Acts::VectorHelpers::phi(track.momentum()));
  }

  std::vector<std::size_t> cleanTracks =
      getCleanedOutTracks(trackScore, counterMaps, measurementsPerTrack);

  ACTS_INFO("Number of tracks: " << tracks.size());
  ACTS_INFO("Number of clean tracks: " << cleanTracks.size());
  ACTS_INFO("Min score: " << m_cfg.minScore);

  std::vector<int> goodTracks;
  std::size_t iTrack = 0;
  for (const auto& track : tracks) {
    if (std::find(cleanTracks.begin(), cleanTracks.end(), iTrack) !=
        cleanTracks.end()) {
      if (trackScore[iTrack] >= m_cfg.minScore) {
        goodTracks.push_back(track.index());
      }
    }
    iTrack++;
  }

  ACTS_INFO("Number of good tracks: " << goodTracks.size());
  return goodTracks;
}

}  // namespace Acts
