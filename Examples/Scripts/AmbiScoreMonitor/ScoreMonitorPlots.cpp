// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
 
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>

// Enum for match types
enum MatchType {
    GOOD,
    FAKE,
    DUPLICATE,
    UNKNOWN,
};

// Helper function to determine match type
MatchType getMatchType(bool isFake, bool isDuplicate, bool isGood) {
    std::cout << "isFake: " << isFake << ", isDuplicate: " << isDuplicate 
              << ", isGood: " << isGood << std::endl;
    if (isGood) return GOOD;
    if (isFake) return FAKE;
    if (isDuplicate) return DUPLICATE;
    return UNKNOWN;
}

// Helper function to get color for match type
int getMatchColor(MatchType type) {
    switch (type) {
        case GOOD: return kGreen;
        case FAKE: return kRed;
        case DUPLICATE: return kOrange;
        default: return kGray;
    }
}

// Helper function to get name for match type
std::string getMatchName(MatchType type) {
    switch (type) {
        case GOOD: return "Good Match";
        case FAKE: return "Fake";
        case DUPLICATE: return "Duplicate";
        default: return "Unknown";
    }
}

// Helper to create and draw scatter plot with match type coloring
TGraph *makeScatterByMatchType(const std::vector<double> &x,
                              const std::vector<double> &y,
                              const std::vector<MatchType> &matchTypes,
                              const std::string &title,
                              int markerStyle = 20)
{
    if (x.empty() || y.empty() || x.size() != y.size() || x.size() != matchTypes.size()) {
        std::cerr << "Error: Vector sizes don't match or are empty for: " << title << std::endl;
        return nullptr;
    }

    // for (size_t i = 0; i < matchTypes.size(); i++) {
    //   std::cout << "MatchType" << i << ": " << getMatchName(matchTypes[i]) << std::endl;
    // }
    
    // Separate points by match type
    std::vector<std::vector<double>> x_by_type(4); // 4 match types
    std::vector<std::vector<double>> y_by_type(4);
    
    for (size_t i = 0; i < x.size(); i++) {
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
            TGraph *g = new TGraph(x_by_type[type].size(), x_by_type[type].data(), y_by_type[type].data());
            g->SetMarkerStyle(markerStyle);
            g->SetMarkerColor(getMatchColor(static_cast<MatchType>(type)));

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

// Function to process tree and extract data
bool processTree(TTree *tree, 
                std::vector<double> &v_pt, 
                std::vector<double> &v_eta, 
                std::vector<double> &v_chi2, 
                std::vector<double> &v_total,
                std::vector<MatchType> &v_matchTypes,
                std::vector<std::vector<double>> &v_hit,
                std::vector<std::vector<double>> &v_hole,
                std::vector<std::string> &v_detectorNames)
{
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

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Only process matched tracks
        if (!isMatched) continue;
        
        v_pt.push_back(pt);
        v_eta.push_back(eta);
        v_chi2.push_back(chi2Score);
        v_total.push_back(totalScore);
        
        MatchType matchType = getMatchType(isFake, isDuplicate, isGood);
        v_matchTypes.push_back(matchType);

        // Process detector scores
        if (detectorHitScore && !detectorHitScore->empty()) {
            if (v_hit.empty()) v_hit.resize(detectorHitScore->size());
            for (size_t j = 0; j < detectorHitScore->size() && j < v_hit.size(); j++) {
                v_hit[j].push_back(detectorHitScore->at(j));
            }
        }
        
        if (detectorHoleScore && !detectorHoleScore->empty()) {
            if (v_hole.empty()) v_hole.resize(detectorHoleScore->size());
            for (size_t j = 0; j < detectorHoleScore->size() && j < v_hole.size(); j++) {
                v_hole[j].push_back(detectorHoleScore->at(j));
            }
        }
        
        if (detectorNames && !detectorNames->empty() && v_detectorNames.empty()) {
            v_detectorNames = *detectorNames;
        }
    }
    
    return true;
}

// Function to create basic scatter plots
void createBasicPlots(const std::vector<double> &v_pt,
                     const std::vector<double> &v_eta,
                     const std::vector<double> &v_chi2,
                     const std::vector<double> &v_total,
                     const std::vector<MatchType> &v_matchTypes)
{
    if (v_pt.empty()) {
        std::cerr << "No matched tracks found for basic plots!" << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "Basic Scatter Plots", 1600, 900);
    c1->Divide(2, 2);

    std::vector<std::vector<double>> x_axis = {v_pt, v_eta};
    std::vector<std::string> xlabels = {"pT", "eta"};

    for (size_t i = 0; i < 2; ++i) {
        c1->cd(i + 1);
        makeScatterByMatchType(x_axis[i], v_chi2, v_matchTypes, xlabels[i] + " vs chi2Score", 21);
        
        c1->cd(3 + i);
        makeScatterByMatchType(x_axis[i], v_total, v_matchTypes, xlabels[i] + " vs totalScore", 21);
    }

    c1->SaveAs("scatter_pt_eta_matched.png");
    std::cout << "Saved: scatter_pt_eta_matched.png" << std::endl;
    delete c1;
}

// Function to create detector-specific plots
void createDetectorPlots(const std::vector<double> &v_pt,
                        const std::vector<double> &v_eta,
                        const std::vector<MatchType> &v_matchTypes,
                        const std::vector<std::vector<double>> &v_scores,
                        const std::vector<std::string> &v_detectorNames,
                        const std::string &scoreType,
                        const std::string &filename)
{
    if (v_scores.empty()) {
        std::cerr << "No " << scoreType << " data available!" << std::endl;
        return;
    }
    
    TCanvas *c = new TCanvas(("c_" + scoreType).c_str(), (scoreType + " Plots").c_str(), 1600, 900);
    c->Divide(v_scores.size(), 2);

    for (size_t i = 0; i < v_scores.size(); ++i) {
        if (v_scores[i].empty()) continue;
        
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
    
    c->SaveAs(filename.c_str());
    std::cout << "Saved: " << filename << std::endl;
    delete c;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.root" << std::endl;
        return 1;
    }

    std::cout << "Opening file: " << argv[1] << std::endl;
    TFile *file = TFile::Open(argv[1], "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 1;
    }

    TTree *tree = dynamic_cast<TTree*>(file->Get("scoremonitor"));
    if (!tree) {
        std::cerr << "Could not find tree scoremonitor" << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    std::cout << "Found tree with " << tree->GetEntries() << " entries" << std::endl;

    // Storage vectors
    std::vector<double> v_pt, v_eta, v_chi2, v_total;
    std::vector<MatchType> v_matchTypes;
    std::vector<std::vector<double>> v_hit, v_hole;
    std::vector<std::string> v_detectorNames;

    // Process the tree
    if (!processTree(tree, v_pt, v_eta, v_chi2, v_total, v_matchTypes, v_hit, v_hole, v_detectorNames)) {
        std::cerr << "Error processing tree!" << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    std::cout << "Processed " << v_pt.size() << " matched tracks" << std::endl;

    // Create plots
    createBasicPlots(v_pt, v_eta, v_chi2, v_total, v_matchTypes);
    createDetectorPlots(v_pt, v_eta, v_matchTypes, v_hit, v_detectorNames, "Hit", "scatter_hit_scores_matched.png");
    createDetectorPlots(v_pt, v_eta, v_matchTypes, v_hole, v_detectorNames, "Hole", "scatter_hole_scores_matched.png");

    // Close the file
    file->Close();
    delete file;
    
    std::cout << "Program completed successfully!" << std::endl;
    return 0;
}