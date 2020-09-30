#include "AnalysisTree/Cuts.hpp"
#include "CbmCuts.h"

#include "src/EntryConfig.hpp"
#include "src/Manager.hpp"
#include "src/Task.hpp"
#include "src/Utils.hpp"

using namespace AnalysisTree;

void BuildEffmap(QA::Task& task);

std::string sim_particles = "SimParticles";
std::string reco_tracks = "LambdaCandidates";
const float y_beam = 1.62179;
const float HugeValue = 9.9e9;

SimpleCut pid_cut({sim_particles, "pid"}, 3122);
Cuts* sim_pid_cut = new Cuts("LambdaSimCut", {pid_cut});

SimpleCut signal_cut({reco_tracks, "is_signal"}, 1, 2);
SimpleCut cosinepos_cut({reco_tracks, "cosinepos"}, 0.9999, HugeValue);
Cuts* selection_cuts = new Cuts("LambdaCandidatesCuts", {
//                                                           nhitspos_cut,
//                                                           nhitsneg_cut,
//                                                           sumnhits_cut,
//                                                           chi2primpos_cut,
//                                                           chi2primneg_cut,
                                                          cosinepos_cut,
//                                                           ldl_cut,
                                                          signal_cut
                                                                          });

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./pfs_qa filelist1 filelist2" << std::endl;
    return -1;
  }

  const std::string filelist1 = argv[1];
  const std::string filelist2 = argv[2];
  
  QA::Manager man({filelist1, filelist2}, {"aTree", "sTree"});
  man.SetOutFileName("effmap.root");
  
  man.AddBranchCut(sim_pid_cut);
  man.AddBranchCut(selection_cuts);
  
  auto* task = new QA::Task;

  BuildEffmap(*task);
  
  man.AddTask(task);

  man.Init();
  man.Run(-1);
  man.Finish();

  return 0;
}

void BuildEffmap(QA::Task& task) {
  
  const int y_nbins = 22;
  const float y_low = y_beam-0.8;
  const float y_up = y_beam+1.4;
  
  const int phi_nbins = 30;
  const float phi_low = -TMath::Pi();
  const float phi_up = TMath::Pi();
  
  const int pT_nbins = 25;
  const float pT_low = 0.;
  const float pT_up = 2.5;
  
  task.AddH2({"y_{LAB}", {sim_particles, "rapidity"}, {y_nbins, y_low, y_up}}, {"p_{T}, GeV/c", {sim_particles, "pT"}, {pT_nbins, pT_low, pT_up}});
  task.AddH2({"#varphi", {sim_particles, "phi"}, {phi_nbins, phi_low, phi_up}}, {"p_{T}, GeV/c", {sim_particles, "pT"}, {pT_nbins, pT_low, pT_up}});
  task.AddH2({"y_{LAB}", {reco_tracks, "rapidity"}, {y_nbins, y_low, y_up}}, {"p_{T}, GeV/c", {reco_tracks, "pT"}, {pT_nbins, pT_low, pT_up}});
  task.AddH2({"#varphi", {reco_tracks, "phi"}, {phi_nbins, phi_low, phi_up}}, {"p_{T}, GeV/c", {reco_tracks, "pT"}, {pT_nbins, pT_low, pT_up}});
}
