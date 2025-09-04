// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <array>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootScoreMonitorWriter
///
/// Write out tracks (i.e. a vector of trackState at the moment) into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to one track for optimum writing speed.
/// The event number is part of the written data.
///
/// A common file can be provided for the writer to attach his TTree, this is
/// done by setting the Config::rootFile pointer to an existing file.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootScoreMonitorWriter final
    : public WriterT<
          std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>> {
 public:
  struct Config {
    /// Input particles collection.
    std::string inputScoreMonitor;

    /// Input particles collection.
    std::string inputParticles;
    
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;

    /// Names of the detectors
    std::vector<std::string> detectorNames = {};
    /// Names of optional scoring functions
    std::vector<std::string> optionalFunctions = {};

    /// output filename.
    std::string filePath = "ScoreMonitor.root";
    /// name of the output tree.
    std::string treeName = "ScoreMonitorTree";
    /// file access mode.
    std::string fileMode = "RECREATE";
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootScoreMonitorWriter(const Config& config, Acts::Logging::Level level);

  ~RootScoreMonitorWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] tracks are what to be written out
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>&
          scoreMonitor) override;

 private:
  /// The config class
  Config m_cfg;

  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  /// Mutex used to protect multi-threaded writes
  std::mutex m_writeMutex;

  /// The output file
  TFile* m_outputFile{nullptr};
  /// The output tree
  TTree* m_outputTree{nullptr};

  /// the event number
  std::uint32_t m_eventNr{0};
  /// transverse momentum
  double m_pT = 0.0;
  /// pseudorapidity
  double m_eta = 0.0;
  /// azimuthal angle
  double m_phi = 0.0;
  /// track index
  std::uint32_t m_index = -1;
  /// total score from ambiguity resolution
  double m_totalScore = 0;
  /// chi2 score from ambiguity resolution
  double m_chi2Score = 0;
  /// pT score from ambiguity resolution
  double m_ptScore = 0;
  /// whether the track is matched to a truth particle
  bool m_isMatched = false;
  /// whether the track is a duplicate
  bool m_isDuplicate = false;
  /// whether the track is a fake
  bool m_isFake = false;
  /// whether the track is a good match
  bool m_isGood = false;

  /// ID of the matched truth particle (0 if not matched)
  std::uint32_t m_matchedParticleId = 0;

  /// scores per detector for hits, holes, outliers, others
  /// (size is number of detectors)
  std::vector<double> m_detectorHitScore;
  std::vector<double> m_detectorHoleScore;
  std::vector<double> m_detectorOutlierScore;
  std::vector<double> m_detectorOtherScore;

  std::vector<double> m_optionalScore;
  std::vector<std::string> m_detectorNamesroot;
  std::vector<std::string> m_optionalFunctionsroot;
};

}  // namespace ActsExamples
