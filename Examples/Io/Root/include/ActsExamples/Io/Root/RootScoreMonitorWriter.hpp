// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"
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
class RootScoreMonitorWriter final : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (fitted) tracks collection
    std::string inputTracks;
    /// Input particles collection.
    std::string inputScoreMonitor;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;

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
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

 private:
  /// The config class
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>>
      m_outputScoreMonitor{this, "InputScoreMonitor"};
  /// Mutex used to protect multi-threaded writes
  std::mutex m_writeMutex;

  /// The output file
  TFile* m_outputFile{nullptr};
  /// The output tree
  TTree* m_outputTree{nullptr};

  /// the event number
  std::uint32_t m_eventNr{0};
  /// the track number
  std::uint32_t m_trackNr{0};

  double m_pT = 0.0;
  double m_eta = 0.0;
  double m_phi = 0.0;
  std::uint32_t m_index = -1;
  double m_totalScore = 0;
  double m_chi2Score = 0;
  double m_ptScore = 0;

  std::vector<double> m_detectorHitScore;
  std::vector<double> m_detectorHoleScore;
  std::vector<double> m_detectorOutlierScore;
  std::vector<double> m_detectorOtherScore;

  std::vector<double> m_optionalScore;
  std::vector<std::string> m_detectorNamesroot;
  std::vector<std::string> m_optionalFunctionsroot;
};

}  // namespace ActsExamples
