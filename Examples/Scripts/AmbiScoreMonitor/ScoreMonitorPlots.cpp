// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TTree.h"

// Enum for match types
enum MatchType {
  GOOD,
  FAKE,
  DUPLICATE,
  UNKNOWN,
};

// Helper function to determine match type
MatchType getMatchType(bool isFake, bool isDuplicate, bool isGood) {
  if (isGood)
    return GOOD;
  if (isFake)
    return FAKE;
  if (isDuplicate)
    return DUPLICATE;
  return UNKNOWN;
}

// Helper function to get color for match type
int getMatchColor(MatchType type) {
  switch (type) {
    case GOOD:
      return kGreen + 2;
    case FAKE:
      return kRed + 1;
    case DUPLICATE:
      return kOrange + 7;
    default:
      return kGray;
  }
}

// Helper function to get name for match type
std::string getMatchName(MatchType type) {
  switch (type) {
    case GOOD:
      return "Good Match";
    case FAKE:
      return "Fake";
    case DUPLICATE:
      return "Duplicate";
    default:
      return "Unknown";
  }
}

// Helper to create and draw scatter plot with match type coloring
TGraph *makeScatterByMatchType(const std::vector<double> &x,
                               const std::vector<double> &y,
                               const std::vector<MatchType> &matchTypes,
                               const std::string &title, int markerStyle = 20) {
  if (x.empty() || y.empty() || x.size() != y.size() ||
      x.size() != matchTypes.size()) {
    std::cout << "size x: " << x.size() << ", size y: " << y.size()
              << ", size matchTypes: " << matchTypes.size() << std::endl;
    std::cerr << "Error: Vector sizes don't match or are empty for: " << title
              << std::endl;
    return nullptr;
  }

  // Separate points by match type
  std::vector<std::vector<double>> x_by_type(4);  // 4 match types
  std::vector<std::vector<double>> y_by_type(4);

  for (std::size_t i = 0; i < x.size(); i++) {
    int type_idx = static_cast<int>(matchTypes[i]);
    if (type_idx >= 0 && type_idx < 5) {
      x_by_type[type_idx].push_back(x[i]);
      y_by_type[type_idx].push_back(y[i]);
    }
  }

  TGraph *first_graph = nullptr;
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

  // Create graphs for each match type
  for (int type = 0; type < 4; type++) {
    if (!x_by_type[type].empty()) {
      TGraph *g = new TGraph(x_by_type[type].size(), x_by_type[type].data(),
                             y_by_type[type].data());
      g->SetMarkerStyle(markerStyle);

      int baseColor = getMatchColor(static_cast<MatchType>(type));
      g->SetMarkerColorAlpha(baseColor, 0.5);
      if (!first_graph) {
        g->SetTitle(title.c_str());
        g->Draw("AP");
        first_graph = g;
      } else {
        g->Draw("P SAME");
      }

      leg->AddEntry(g, getMatchName(static_cast<MatchType>(type)).c_str(), "p");
    }
  }

  if (first_graph) {
    leg->Draw();
  } else {
    delete leg;
  }

  return first_graph;
}
void processDetectorScores(std::vector<double> *detectorScores,
                           std::vector<std::vector<double>> &targetVector) {
  if (detectorScores) {
    if (targetVector.empty()) {
      targetVector.resize(detectorScores->size());
    }
    for (std::size_t j = 0;
         j < detectorScores->size() && j < targetVector.size(); j++) {
      targetVector[j].push_back(detectorScores->at(j));
    }
  }
}

// Function to process tree and extract data
bool processTree(TTree *tree, std::vector<double> &v_pt,
                 std::vector<double> &v_eta, std::vector<double> &v_chi2,
                 std::vector<double> &v_total,
                 std::vector<MatchType> &v_matchTypes,
                 std::vector<std::vector<double>> &v_hit,
                 std::vector<std::vector<double>> &v_hole,
                 std::vector<std::vector<double>> &v_outlier,
                 std::vector<std::vector<double>> &v_other,
                 std::vector<std::vector<double>> &v_optional,
                 std::vector<std::string> &v_detectorNames) {
  // Variables
  double pt = 0, eta = 0, phi = 0;
  double ptScore = 0, chi2Score = 0, totalScore = 0;
  bool isMatched = false;
  bool isFake = false, isDuplicate = false, isGood = false;
  std::vector<double> *detectorHitScore = nullptr;
  std::vector<double> *detectorHoleScore = nullptr;
  std::vector<double> *detectorOutlierScore = nullptr;
  std::vector<double> *detectorOtherScore = nullptr;
  std::vector<double> *optionalScore = nullptr;
  std::vector<std::string> *detectorNames = nullptr;

  // Set branch addresses
  tree->SetBranchAddress("pT", &pt);
  tree->SetBranchAddress("eta", &eta);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("isFake", &isFake);
  tree->SetBranchAddress("isDuplicate", &isDuplicate);
  tree->SetBranchAddress("isGood", &isGood);
  tree->SetBranchAddress("ptScore", &ptScore);
  tree->SetBranchAddress("chi2Score", &chi2Score);
  tree->SetBranchAddress("totalScore", &totalScore);
  tree->SetBranchAddress("isMatched", &isMatched);

  // Detector branches
  if (tree->GetBranch("detectorHitScore")) {
    tree->SetBranchAddress("detectorHitScore", &detectorHitScore);
  }
  if (tree->GetBranch("detectorHoleScore")) {
    tree->SetBranchAddress("detectorHoleScore", &detectorHoleScore);
  }
  if (tree->GetBranch("detectorNames")) {
    tree->SetBranchAddress("detectorNames", &detectorNames);
  }
  if (tree->GetBranch("detectorOutlierScore")) {
    tree->SetBranchAddress("detectorOutlierScore", &detectorOutlierScore);
  }
  if (tree->GetBranch("detectorOtherScore")) {
    tree->SetBranchAddress("detectorOtherScore", &detectorOtherScore);
  }
  if (tree->GetBranch("optionalScore")) {
    tree->SetBranchAddress("optionalScore", &optionalScore);
  }

  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);

    // Only process matched tracks
    if (!isMatched)
      continue;
    if (detectorHitScore->empty())
      continue;
    if (detectorHoleScore->empty())
      continue;

    v_pt.push_back(pt);
    v_eta.push_back(eta);
    v_chi2.push_back(chi2Score);
    v_total.push_back(totalScore);

    MatchType matchType = getMatchType(isFake, isDuplicate, isGood);
    v_matchTypes.push_back(matchType);

    processDetectorScores(detectorHitScore, v_hit);
    processDetectorScores(detectorHoleScore, v_hole);
    processDetectorScores(detectorOutlierScore, v_outlier);
    processDetectorScores(detectorOtherScore, v_other);
    processDetectorScores(optionalScore, v_optional);

    if (detectorNames && !detectorNames->empty() && v_detectorNames.empty()) {
      v_detectorNames = *detectorNames;
    }
  }

  return true;
}

// Function to create stacked histogram of total scores
void createStackedTotalScoreHistogram(
    const std::vector<double> &v_total,
    const std::vector<MatchType> &v_matchTypes, const std::string &outputDir,
    const std::string &outputExt) {
  if (v_total.empty()) {
    std::cerr << "No data found for stacked histogram!" << std::endl;
    return;
  }

  // Create output path
  std::string outputPath = outputDir + "stacked_total_score" + outputExt;

  // Find min and max of total scores for binning
  double min_score = *std::min_element(v_total.begin(), v_total.end());
  double max_score = *std::max_element(v_total.begin(), v_total.end());

  // Create histogram with reasonable binning
  int n_bins = 50;
  THStack *stack =
      new THStack("hstack", "Total Score Distribution;Total Score;Counts");

  // Create histograms for each category
  TH1F *h_good = new TH1F("h_good", "Good", n_bins, min_score, max_score);
  TH1F *h_duplicate =
      new TH1F("h_duplicate", "Duplicate", n_bins, min_score, max_score);
  TH1F *h_fake = new TH1F("h_fake", "Fake", n_bins, min_score, max_score);
  TH1F *h_unknown =
      new TH1F("h_unknown", "Unknown", n_bins, min_score, max_score);

  // Set colors using your helper function
  h_good->SetFillColorAlpha(getMatchColor(GOOD), 0.8);
  h_duplicate->SetFillColorAlpha(getMatchColor(DUPLICATE), 0.8);
  h_fake->SetFillColorAlpha(getMatchColor(FAKE), 0.8);
  h_unknown->SetFillColorAlpha(getMatchColor(UNKNOWN), 0.8);

  // Fill histograms based on match type
  for (std::size_t i = 0; i < v_total.size(); ++i) {
    switch (v_matchTypes[i]) {
      case GOOD:
        h_good->Fill(v_total[i]);
        break;
      case DUPLICATE:
        h_duplicate->Fill(v_total[i]);
        break;
      case FAKE:
        h_fake->Fill(v_total[i]);
        break;
      case UNKNOWN:
      default:
        h_unknown->Fill(v_total[i]);
        break;
    }
  }

  // Add to stack in order (background first)
  stack->Add(h_unknown);
  stack->Add(h_fake);
  stack->Add(h_duplicate);
  stack->Add(h_good);

  // Create canvas and draw
  TCanvas *c = new TCanvas("c_stack", "Stacked Total Score", 800, 600);
  stack->Draw("HIST");

  // Add legend
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(h_good, "Good", "f");
  leg->AddEntry(h_duplicate, "Duplicate", "f");
  leg->AddEntry(h_fake, "Fake", "f");
  if (h_unknown->GetEntries() > 0) {
    leg->AddEntry(h_unknown, "Unknown", "f");
  }
  leg->Draw();

  // Add statistics
  TPaveText *stats = new TPaveText(0.15, 0.7, 0.4, 0.9, "NDC");
  stats->SetFillColor(0);
  stats->SetBorderSize(1);
  stats->AddText(Form("Good: %d", static_cast<int>(h_good->GetEntries())));
  stats->AddText(
      Form("Duplicate: %d", static_cast<int>(h_duplicate->GetEntries())));
  stats->AddText(Form("Fake: %d", static_cast<int>(h_fake->GetEntries())));
  stats->AddText(
      Form("Unknown: %d", static_cast<int>(h_unknown->GetEntries())));
  stats->Draw();

  c->SaveAs(outputPath.c_str());
  std::cout << "Saved: " << outputPath << std::endl;

  // Cleanup
  delete c;
  delete stack;
  delete h_good;
  delete h_duplicate;
  delete h_fake;
  delete h_unknown;
  delete leg;
  delete stats;
}

// Function to create basic scatter plots
void createBasicPlots(const std::vector<double> &v_pt,
                      const std::vector<double> &v_eta,
                      const std::vector<double> &v_chi2,
                      const std::vector<double> &v_total,
                      const std::vector<MatchType> &v_matchTypes,
                      const std::string &outputDir,
                      const std::string &outputExt) {
  if (v_pt.empty()) {
    std::cerr << "No matched tracks found for basic plots!" << std::endl;
    return;
  }

  // Create output path
  std::string outputPath = outputDir + "scatter_pt_eta_matched" + outputExt;

  TCanvas *c1 = new TCanvas("c1", "Basic Scatter Plots", 1600, 900);
  c1->Divide(2, 2);

  std::vector<std::vector<double>> x_axis = {v_pt, v_eta};
  std::vector<std::string> xlabels = {"pT", "eta"};

  for (std::size_t i = 0; i < 2; ++i) {
    c1->cd(i + 1);
    makeScatterByMatchType(x_axis[i], v_chi2, v_matchTypes,
                           xlabels[i] + " vs chi2Score", 21);

    c1->cd(3 + i);
    makeScatterByMatchType(x_axis[i], v_total, v_matchTypes,
                           xlabels[i] + " vs totalScore", 21);
  }

  c1->SaveAs(outputPath.c_str());
  std::cout << "Saved: " << outputPath << std::endl;
  delete c1;
}

// Function to create detector-specific plots
void createDetectorPlots(const std::vector<double> &v_pt,
                         const std::vector<double> &v_eta,
                         const std::vector<MatchType> &v_matchTypes,
                         const std::vector<std::vector<double>> &v_scores,
                         const std::vector<std::string> &v_detectorNames,
                         const std::string &scoreType,
                         const std::string &outputDir,
                         const std::string &outputExt) {
  if (v_scores.empty()) {
    std::cerr << "No " << scoreType << " data available!" << std::endl;
    return;
  }

  std::string outputPath =
      outputDir + "scatter_" + scoreType + "_scores_matched" + outputExt;

  TCanvas *c = new TCanvas(("c_" + scoreType).c_str(),
                           (scoreType + " Plots").c_str(), 1600, 900);
  c->Divide(v_scores.size(), 2);

  for (std::size_t i = 0; i < v_scores.size(); ++i) {
    if (v_scores[i].empty())
      continue;

    std::string detectorLabel;
    if (!v_detectorNames.empty() && i < v_detectorNames.size()) {
      detectorLabel = v_detectorNames[i];
    } else {
      std::stringstream ss;
      ss << "Detector_" << i;
      detectorLabel = ss.str();
    }

    c->cd(i + 1);
    makeScatterByMatchType(v_pt, v_scores[i], v_matchTypes,
                           "pT vs " + scoreType + " " + detectorLabel, 20);

    c->cd(v_scores.size() + i + 1);
    makeScatterByMatchType(v_eta, v_scores[i], v_matchTypes,
                           "eta vs " + scoreType + " " + detectorLabel, 20);
  }
  c->SaveAs(outputPath.c_str());

  std::cout << "Saved: " << outputPath << std::endl;
  delete c;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " input.root [output_directory] [output_extension]"
              << std::endl;
    std::cerr
        << "Default: output_directory = ./output/, output_extension = .pdf"
        << std::endl;
    return 1;
  }

  // Set default output parameters
  std::string outputDir = "./output/";
  std::string outputExt = ".pdf";

  // Parse command line arguments
  if (argc >= 3) {
    outputDir = argv[2];
    // Ensure output directory ends with a slash
    if (outputDir.back() != '/') {
      outputDir += "/";
    }
  }
  if (argc >= 4) {
    outputExt = argv[3];
    // Ensure extension starts with a dot
    if (outputExt[0] != '.') {
      outputExt = "." + outputExt;
    }
  }

  std::cout << "Opening file: " << argv[1] << std::endl;
  std::cout << "Output directory: " << outputDir << std::endl;
  std::cout << "Output extension: " << outputExt << std::endl;

  TFile *file = TFile::Open(argv[1], "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << argv[1] << std::endl;
    return 1;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("scoremonitor"));
  if (!tree) {
    std::cerr << "Could not find tree scoremonitor" << std::endl;
    file->Close();
    delete file;
    return 1;
  }

  std::cout << "Found tree with " << tree->GetEntries() << " entries"
            << std::endl;

  // Storage vectors
  std::vector<double> v_pt, v_eta, v_chi2, v_total;
  std::vector<MatchType> v_matchTypes;
  std::vector<std::vector<double>> v_hit, v_hole, v_outlier, v_other,
      v_optional;
  std::vector<std::string> v_detectorNames;

  // Process the tree
  if (!processTree(tree, v_pt, v_eta, v_chi2, v_total, v_matchTypes, v_hit,
                   v_hole, v_outlier, v_other, v_optional, v_detectorNames)) {
    std::cerr << "Error processing tree!" << std::endl;
    file->Close();
    delete file;
    return 1;
  }

  std::cout << "Processed " << v_pt.size() << " matched tracks" << std::endl;

  /// Create stacked total score histogram
  createStackedTotalScoreHistogram(v_total, v_matchTypes, outputDir, outputExt);

  // Create plots with proper output paths
  createBasicPlots(v_pt, v_eta, v_chi2, v_total, v_matchTypes, outputDir,
                   outputExt);
  createDetectorPlots(v_pt, v_eta, v_matchTypes, v_hit, v_detectorNames, "Hit",
                      outputDir, outputExt);
  createDetectorPlots(v_pt, v_eta, v_matchTypes, v_hole, v_detectorNames,
                      "Hole", outputDir, outputExt);

  // Close the file
  file->Close();
  delete file;

  std::cout << "Program completed successfully!" << std::endl;
  return 0;
}