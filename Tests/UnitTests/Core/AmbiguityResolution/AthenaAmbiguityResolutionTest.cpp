// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/EventData/detail/TestTrackState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <map>

using Acts::MultiTrajectoryTraits::IndexType;

namespace Acts::detail::Test {
VectorMultiTrajectory create() {
  return {};
}
ConstVectorMultiTrajectory createConst() {
  return {};
}

TrackContainer<VectorTrackContainer, VectorMultiTrajectory, RefHolder>
createfaketracks(std::size_t nTracks) {
  using trajectory_t = VectorMultiTrajectory;

  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  std::vector<trajectory_t> trajectories;
  std::vector<TestTrackState> trackStates;
  auto mkts = [&](auto prev) {
    std::default_random_engine rng;
    VectorMultiTrajectory& traj = mtj;
    if constexpr (std::is_same_v<decltype(prev), IndexType>) {
      auto ts =
          traj.getTrackState(traj.addTrackState(TrackStatePropMask::All, prev));
      TestTrackState pc(rng, 2u);
      fillTrackState<VectorMultiTrajectory>(pc, TrackStatePropMask::All, ts);
      return ts;
    } else {
      auto ts = traj.getTrackState(
          traj.addTrackState(TrackStatePropMask::All, prev.index()));
      TestTrackState pc(rng, 2u);
      fillTrackState<VectorMultiTrajectory>(pc, TrackStatePropMask::All, ts);
      return ts;
    }
  };

  for (unsigned int i = 0; i < nTracks; i++) {
    auto ts1 = mkts(MultiTrajectoryTraits::kInvalid);
    auto ts2 = mkts(ts1);
    auto ts3 = mkts(ts2);
    auto ts4 = mkts(ts3);
    auto ts5 = mkts(ts4);
    auto ts6 = mkts(ts5);
    auto ts7 = mkts(ts6);
    auto ts8 = mkts(ts7);
    auto ts9 = mkts(ts8);
    auto ts10 = mkts(ts9);

    auto t = tc.makeTrack();
    t.tipIndex() = ts10.index();
  }
  // BOOST_CHECK_EQUAL(tc.size(), nTracks);

  unsigned int i = 0;
  for (auto track : tc) {
    // BOOST_CHECK_EQUAL(i, track.tipIndex());
    track.parameters().setRandom();
    i++;
  }

  // BOOST_CHECK_EQUAL(std::distance(tc.begin(), tc.end()), tc.size());

  return tc;
}
}  // namespace Acts::detail::Test

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(AthenaAmbiguityResolutionTest)

// Test fixture for AthenaAmbiguityResolution
struct Fixture {
  AthenaAmbiguityResolution::Config config;

  using Counter = AthenaAmbiguityResolution::Counter;
  using DetectorConfig = AthenaAmbiguityResolution::DetectorConfig;

  Fixture() {
    // Set up any resources used by the tests
    config.volumeMap = {{8, 0},  {9, 0},  {10, 0}, {13, 0}, {14, 0},
                        {15, 0}, {16, 0}, {17, 0}, {18, 0}, {19, 0},
                        {20, 0}, {22, 1}, {23, 1}, {24, 1}};

    auto tempDetector = DetectorConfig();
    std::map<std::size_t, DetectorConfig> detectorMap = {{0, tempDetector},
                                                         {1, tempDetector}};
    config.detectorMap = detectorMap;

    config.minScore = 0;
    config.minScoreSharedTracks = 0;
    config.maxShared = 5;
    config.maxSharedTracksPerMeasurement = 10;
    config.phiMax = 3.14;
    config.phiMin = -3.14;
    config.etaMax = 2.7;
    config.etaMin = -2.7;
    config.pTMax = 1400;
    config.pTMin = 0.5;
    config.useAmbigFcn = false;
  }

  ~Fixture() = default;
};

// Helper function to create a sample input for getCleanedOutTracks
std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
createSampleInput() {
  std::vector<std::pair<std::size_t, std::vector<std::size_t>>> trackVolumes = {
      {0, {19, 18, 18, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
      {1, {19, 18, 18, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
      {2, {13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8}},
      {3, {13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8}},
      {4, {19, 18, 18, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10}}};

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
      measurementsPerTrack;
  // Add sample measurements for each track

  for (const auto& trackVolume : trackVolumes) {
    std::vector<std::tuple<std::size_t, std::size_t, bool>> measurements;
    for (std::size_t i = 0; i < trackVolume.second.size(); ++i) {
      measurements.push_back(std::make_tuple(trackVolume.second[i], i, false));
    }
    measurementsPerTrack.push_back(measurements);
  }

  return measurementsPerTrack;
}

BOOST_FIXTURE_TEST_CASE(simpleScoreTest, Fixture) {
  Fixture fixture;
  // Create an instance of AthenaAmbiguityResolution

  AthenaAmbiguityResolution tester(fixture.config);

  TrackContainer<VectorTrackContainer, VectorMultiTrajectory, detail::RefHolder>
      tracks = detail::Test::createfaketracks(5);

  // Assert the expected results
  BOOST_CHECK_EQUAL(tracks.size(), 5);
  // BOOST_CHECK_EQUAL(score[0], 0);
  // BOOST_CHECK_EQUAL(score[1], 0);
  // BOOST_CHECK_EQUAL(score[2], 0);
}

BOOST_FIXTURE_TEST_CASE(GetCleanedOutTracksTest, Fixture) {
  Fixture fixture;
  // Create an instance of AthenaAmbiguityResolution
  AthenaAmbiguityResolution tester(fixture.config);

  // Create sample input
  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
      measurementsPerTrack = createSampleInput();

  std::vector<double> TrackSore;
  for (std::size_t i = 0; i < measurementsPerTrack.size(); i++) {
    TrackSore.push_back(60 + 40 * i);
  }

  std::vector<std::map<std::size_t, Counter>> CounterMaps = {
      {{0, {0, 14, 0, 0, 0}}, {1, {0, 0, 0, 0, 0}}},
      {{0, {0, 15, 0, 0, 0}}, {1, {0, 0, 0, 0, 0}}},
      {{0, {0, 17, 0, 0, 0}}, {1, {0, 0, 0, 0, 0}}},
      {{0, {0, 18, 0, 0, 0}}, {1, {0, 0, 0, 0, 0}}},
      {{0, {0, 14, 0, 0, 0}}, {1, {0, 0, 0, 0, 0}}}};

  // Call the function under testBOOST_FIXTURE_TEST_CASE
  std::vector<std::size_t> cleanTracks =
      tester.getCleanedOutTracks(TrackSore, CounterMaps, measurementsPerTrack);

  // Assert the expected results
  BOOST_CHECK_EQUAL(cleanTracks.size(), 5);
  BOOST_CHECK_EQUAL(cleanTracks[0], 0);
  BOOST_CHECK_EQUAL(cleanTracks[1], 1);
  BOOST_CHECK_EQUAL(cleanTracks[2], 2);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
