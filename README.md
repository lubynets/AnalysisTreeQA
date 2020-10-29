# AnalysisTreeQA
[![Build Status](https://travis-ci.com/HeavyIonAnalysis/AnalysisTreeQA.svg?branch=master)](https://travis-ci.com/HeavyIonAnalysis/AnalysisTreeQA)

## Installation

    cd AnalysisTreeQA
    mkdir build install
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install ../
    make -j install
    
## Usage

Simple program: 

    // Manager instance: filelist1, filelist2, ...; treename1, treename2, ...
    QA::Manager man({filelist1, filelist2}, {"aTree", "bTree"}); 
    man.SetOutFileName("test_qa.root");

    auto* task = QA::Task();

    // 1D histo: {branch_name, field_name}, {nbins, min, max}
    task->AddH1({"N_{hits}", {"VtxTracks", "nhits"}, {15, 0, 15}}); 

    // AnalysisTree::Variable in case of more complicated plot
    Variable chi2_over_ndf("chi2/ndf", {{"VtxTracks", "chi2"}, {"VtxTracks", "ndf"}}, []( std::vector<double>& var ) { return var.at(0)/var.at(1); });
    task->AddH1({"#chi^{2}/NDF", chi2_over_ndf, {100, 0, 100}}); 

    // 2D histo:
    task->AddH2({"DCA_{x}", {"VtxTracks", "dcax"}, {100, -1, 1}}, {"DCA_{y}", {"VtxTracks", "dcay"}, {100, -1, 1}}); 

    // Histogram with additional cut:
    Cuts* mc_protons = new Cuts("McProtons", {{"SimTracks","pid"}, 2212});
    task->AddH1({"N_{hits}", {"VtxTracks", "nhits"}, {15, 0, 15}}, mc_protons); 

    man.Run();

## Known problems

 - Matching between 3 branches is not available (i. e. TofHits.m2 vs. VtxTracks.p with cut on MC-true protons) 
