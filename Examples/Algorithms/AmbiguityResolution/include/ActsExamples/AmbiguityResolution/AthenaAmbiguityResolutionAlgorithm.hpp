// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <string>

namespace ActsExamples {

/// Evicts tracks that seem to be duplicated and fake.
///
/// The implementation works as follows:
///  1) Cluster together nearby tracks using shared hits
///  2) For each track use a neural network to compute a score
///  3) In each cluster keep the track with the highest score
class AthenaAmbiguityResolutionAlgorithm final :  public IAlgorithm {
 public:
  /// Configuration for the ambiguity resolution algorithm.

  struct Config {
    /// Input track collection.
    std::string inputTracks;
    /// Output track collection.
    std::string outputTracks;

    std::map<unsigned int,Acts::AthenaAmbiguityResolution::VolumeConfig> volumeMap = {

      {8,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 1
      {9,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 2 (barrel)
      {10,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 0}}, // inner pixel 3

      {13,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 1 
      {14,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 2
      {15,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 3
      {16,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 4 (barrel)
      {18,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 5
      {19,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 6
      {20,{20, -10, 2, 0, 0, 10, 10, 1000, 1000,false, 1}}, // outer pixel 7

      {22,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 1
      {23,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 2 (barrel)
      {24,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // strip 3

      {25,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // HGTD 1
      {2,{10, -5, 2, 0, 0, 10, 10, 1000, 1000,false, 2}}, // HGTD 2
    };
  };


  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  AthenaAmbiguityResolutionAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Acts::AthenaAmbiguityResolution m_core;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
