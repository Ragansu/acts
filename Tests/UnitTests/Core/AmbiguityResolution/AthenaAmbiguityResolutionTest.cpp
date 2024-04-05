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
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFinding/TrackSelector.hpp"
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

  auto makeTrack = [&]() {
    auto track = tc.makeTrack();

    using namespace Acts::UnitLiterals;
    track.parameters() << 0, 0, M_PI / 2, M_PI / 2, 1 / 1_GeV, 0;
    auto perigee = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
    track.setReferenceSurface(perigee);
    return track;
  };

  auto makeSurface = [](GeometryIdentifier id) {
    auto srf =
        Surface::makeShared<PlaneSurface>(Vector3::Zero(), Vector3::UnitZ());

    srf->assignGeometryId(id);
    return srf;
  };

  auto addTrackState = [](auto& track, const auto& surface,
                          TrackStateFlag flag) {
    auto ts = track.appendTrackState();
    ts.setReferenceSurface(surface);
    ts.typeFlags().set(flag);
    return ts;
  };

  auto addMeasurement = [&](auto& track, const auto& surface) {
    return addTrackState(track, surface, TrackStateFlag::MeasurementFlag);
  };

  auto addMaterial = [&](auto& track, const auto& surface) {
    return addTrackState(track, surface, TrackStateFlag::MaterialFlag);
  };

  auto vol7_lay3_sen2 = makeSurface(
      GeometryIdentifier{}.setVolume(7).setLayer(3).setSensitive(2));
  auto vol7_lay4 = makeSurface(GeometryIdentifier{}.setVolume(7).setLayer(4));
  auto vol7_lay3_sen8 = makeSurface(
      GeometryIdentifier{}.setVolume(7).setLayer(3).setSensitive(8));
  auto vol7_lay5_sen11 = makeSurface(
      GeometryIdentifier{}.setVolume(7).setLayer(5).setSensitive(11));
  auto vol7_lay5_sen12 = makeSurface(
      GeometryIdentifier{}.setVolume(7).setLayer(5).setSensitive(12));
  auto vol7_lay6_sen3 = makeSurface(
      GeometryIdentifier{}.setVolume(7).setLayer(6).setSensitive(3));

  auto vol8_lay8_sen1 = makeSurface(
      GeometryIdentifier{}.setVolume(8).setLayer(8).setSensitive(1));
  auto vol8_lay8_sen2 = makeSurface(
      GeometryIdentifier{}.setVolume(8).setLayer(8).setSensitive(2));
  auto vol8_lay9_sen1 = makeSurface(
      GeometryIdentifier{}.setVolume(8).setLayer(9).setSensitive(1));

  for (unsigned int i = 0; i < nTracks; i++) {
    TrackSelector::Config cfgVol7;
    cfgVol7.measurementCounter.addCounter({GeometryIdentifier{}.setVolume(7)},
                                          3);
    TrackSelector selectorVol7{cfgVol7};

    auto trackVol7 = makeTrack();

    BOOST_CHECK(!selectorVol7.isValidTrack(trackVol7));

    // 1 hit in vol7
    addMeasurement(trackVol7, vol7_lay3_sen2);
    addMaterial(trackVol7, vol7_lay4);

    BOOST_CHECK(!selectorVol7.isValidTrack(trackVol7));
    addMeasurement(trackVol7, vol7_lay5_sen11);
    BOOST_CHECK(!selectorVol7.isValidTrack(trackVol7));

    // Now we should have enough hits
    addMeasurement(trackVol7, vol7_lay6_sen3);
    BOOST_CHECK(selectorVol7.isValidTrack(trackVol7));
    
  }
  BOOST_CHECK_EQUAL(tc.size(), nTracks);

  BOOST_CHECK_EQUAL(std::distance(tc.begin(), tc.end()), tc.size());

  return tc;
}

// TrackContainer<VectorTrackContainer, VectorMultiTrajectory, RefHolder>
// createfaketracks_copy(std::size_t nTracks) {
//   // construct initial parameters
//   // create common covariance matrix from reasonable standard deviations
//   Acts::BoundVector stddev;
//   stddev[Acts::eBoundLoc0] = 100_um;
//   stddev[Acts::eBoundLoc1] = 100_um;
//   stddev[Acts::eBoundTime] = 25_ns;
//   stddev[Acts::eBoundPhi] = 2_degree;
//   stddev[Acts::eBoundTheta] = 2_degree;
//   stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
//   Acts::BoundSquareMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
//   // all tracks close to the transverse plane along the x axis w/ small
//   // variations in position, direction.
//   Acts::Vector4 mStartPos0(-3_m, 0.0, 0.0, 1_ns);
//   Acts::Vector4 mStartPos1(-3_m, -15_mm, -15_mm, 2_ns);
//   Acts::Vector4 mStartPos2(-3_m, 15_mm, 15_mm, -1_ns);
//   startParameters = {
//       {mStartPos0, 0_degree, 90_degree, 1_e / 1_GeV, cov, pion},
//       {mStartPos1, -1_degree, 91_degree, 1_e / 1_GeV, cov, pion},
//       {mStartPos2, 1_degree, 89_degree, -1_e / 1_GeV, cov, pion},
//   };
//   Acts::Vector4 mEndPos0(3_m, 0.0, 0.0, 1_ns);
//   Acts::Vector4 mEndPos1(3_m, -100_mm, -100_mm, 2_ns);
//   Acts::Vector4 mEndPos2(3_m, 100_mm, 100_mm, -1_ns);
//   endParameters = {
//       {mEndPos0, 0_degree, 90_degree, 1_e / 1_GeV, cov * 100, pion},
//       {mEndPos1, -1_degree, 91_degree, 1_e / 1_GeV, cov * 100, pion},
//       {mEndPos2, 1_degree, 89_degree, -1_e / 1_GeV, cov * 100, pion},
//   };

//   // create some measurements
//   auto measPropagator = makeStraightPropagator(detector.geometry);
//   std::default_random_engine rng(421235);
//   for (std::size_t trackId = 0u; trackId < startParameters.size(); ++trackId)
//   {
//     auto measurements = createMeasurements(measPropagator, geoCtx, magCtx,
//                                            startParameters[trackId],
//                                            detector.resolutions, rng,
//                                            trackId);
//     for (auto& sl : measurements.sourceLinks) {
//       sourceLinks.emplace(sl.m_geometryId, std::move(sl));
//     }
//   }
// }

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
    config.useAmbiguityFunction = false;
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

  // Call the function under test
  std::vector<std::map<std::size_t, Counter>> counterMaps;
  std::vector<double> score = tester.simpleScore(tracks, counterMaps);

  BOOST_CHECK_EQUAL(score[0], 0);
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
