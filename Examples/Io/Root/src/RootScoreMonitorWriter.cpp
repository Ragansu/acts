// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootScoreMonitorWriter.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ios>
#include <ostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

ActsExamples::RootScoreMonitorWriter::RootScoreMonitorWriter(
    const ActsExamples::RootScoreMonitorWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputScoreMonitor, "RootScoreMonitorWriter", level),
      m_cfg(config) {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // open root file and create the tree
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);

  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("pT", &m_pT);
  m_outputTree->Branch("eta", &m_eta);
  m_outputTree->Branch("phi", &m_phi);
  m_outputTree->Branch("index", &m_index);
  m_outputTree->Branch("ptScore", &m_ptScore);
  m_outputTree->Branch("chi2Score", &m_chi2Score);
  m_outputTree->Branch("totalScore", &m_totalScore);
  m_outputTree->Branch("isMatched", &m_isMatched);
  m_outputTree->Branch("matchedParticleId", &m_matchedParticleId);

  m_outputTree->Branch("detectorHitScore", &m_detectorHitScore);
  m_outputTree->Branch("detectorHoleScore", &m_detectorHoleScore);
  m_outputTree->Branch("detectorOutlierScore", &m_detectorOutlierScore);
  m_outputTree->Branch("detectorOtherScore", &m_detectorOtherScore);
  m_outputTree->Branch("optionalScore", &m_optionalScore);
  m_outputTree->Branch("detectorNames", &m_detectorNamesroot);
  m_outputTree->Branch("optionalFunctions", &m_optionalFunctionsroot);
}

ActsExamples::RootScoreMonitorWriter::~RootScoreMonitorWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootScoreMonitorWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_VERBOSE("Wrote scoreMonitors to tree '" << m_cfg.treeName << "' in '"
                                               << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootScoreMonitorWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>&
        scoreMonitor) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  const auto& particles = m_inputParticles(ctx);

  // Get the event number
  m_eventNr = ctx.eventNumber;
  m_detectorNamesroot = m_cfg.detectorNames;
  m_optionalFunctionsroot = m_cfg.optionalFunctions;

  // Fill the tree
  for (const auto& monitor : scoreMonitor) {
    m_pT = monitor.pT;
    m_eta = monitor.eta;
    m_phi = monitor.phi;
    m_index = monitor.index;

    auto match = trackParticleMatching.find(m_index);
    if (match != trackParticleMatching.end() &&
        match->second.particle.has_value()) {
      if (match != trackParticleMatching.end() &&
          match->second.particle.has_value()) {
        // Get the barcode of the majority truth particle
        auto barcode = match->second.particle.value();
        // Find the truth particle via the barcode
        auto ip = particles.find(barcode);
        if (ip != particles.end()) {
          ACTS_VERBOSE("Find the truth particle with barcode "
                       << barcode << "=" << barcode.value());
        } else {
          ACTS_DEBUG("Truth particle with barcode "
                     << barcode << "=" << barcode.value() << " not found!");
        }
      }
      m_isMatched = true;
      m_matchedParticleId = match->second.particle.value().value();
    } else {
      m_isMatched = false;
      m_matchedParticleId = 0;
    }

    m_ptScore = monitor.ptScore;
    m_chi2Score = monitor.chi2Score;
    m_totalScore = monitor.totalScore;

    m_detectorHitScore = monitor.detectorHitScore;
    m_detectorHoleScore = monitor.detectorHoleScore;
    m_detectorOutlierScore = monitor.detectorOutlierScore;
    m_detectorOtherScore = monitor.detectorOtherScore;
    m_optionalScore = monitor.optionalScore;

    m_outputTree->Fill();
  }
  return ActsExamples::ProcessCode::SUCCESS;
}
