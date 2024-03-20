#include <boost/test/unit_test.hpp>
#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include <map>

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(AthenaAmbiguityResolutionTest)

// Test fixture for AthenaAmbiguityResolution
struct Fixture {

    AthenaAmbiguityResolution::Config config;
    
    Fixture() {
        // Set up any resources used by the tests
        config.volumeMap = {

            {8,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 1
            {9,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 2 (barrel)
            {10,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 3

        };
        config.minScore = 0;
        config.minScoreSharedTracks = 0;
        config.maxShared = 5;
        config.maxSharedTracksPerMeasurement = 10;
        config.phiMax = 3.14;
        config.phiMin = -3.14;
        config.etaMax = 2.7;
        config.etaMin = -2.7;     
   
    }

    ~Fixture() {
        // Clean up any resources used by the tests
    }

    // Helper function to create a sample input for getCleanedOutTracks
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> createSampleInput() {
        std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> measurementsPerTrack;
        // Add sample measurements for each track
        measurementsPerTrack.push_back({{0, 0, false}, {1, 1, false}, {2, 2, false}});
        measurementsPerTrack.push_back({{0, 3, false}, {4, 4, false}, {5, 5, false}});
        measurementsPerTrack.push_back({{6, 6, false}, {7, 7, false}, {8, 8, false}});
        return measurementsPerTrack;
    }
};
BOOST_FIXTURE_TEST_CASE(simpleScoreTest, Fixture) {

    Fixture fixture;
    // Create an instance of AthenaAmbiguityResolution
    AthenaAmbiguityResolution tester = AthenaAmbiguityResolution(fixture.config);

    // Create sample input
    TrackContainer tracks{VectorTrackContainer{},
                            VectorMultiTrajectory{}};
    // Call the function under test
    std::vector<std::map<std::size_t, AthenaAmbiguityResolution::Counter>> counterMaps;

    std::vector<int> score = tester.simpleScore(tracks,counterMaps);

    // Assert the expected results
    BOOST_CHECK_EQUAL(score.size(), 3);
    BOOST_CHECK_EQUAL(score[0], 0);
    BOOST_CHECK_EQUAL(score[1], 0);
    BOOST_CHECK_EQUAL(score[2], 0);
}


BOOST_FIXTURE_TEST_CASE(GetCleanedOutTracksTest, Fixture) {
    Fixture fixture;
    // Create an instance of AthenaAmbiguityResolution
    AthenaAmbiguityResolution tester(fixture.config);

    TrackContainer tracks{VectorTrackContainer{},
                        VectorMultiTrajectory{}};
    // Create sample input
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>> 
    measurementsPerTrack = tester.computeInitialState(tracks, &sourceLinkHash,
                              &sourceLinkEquality);

    // Call the function under testBOOST_FIXTURE_TEST_CASE
    std::vector<std::size_t> cleanTracks = tester.getCleanedOutTracks({}, {}, measurementsPerTrack);

    // Assert the expected results
    BOOST_CHECK_EQUAL(cleanTracks.size(), 3);
    BOOST_CHECK_EQUAL(cleanTracks[0], 0);
    BOOST_CHECK_EQUAL(cleanTracks[1], 1);
    BOOST_CHECK_EQUAL(cleanTracks[2], 2);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
