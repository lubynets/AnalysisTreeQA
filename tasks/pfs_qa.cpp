#include "AnalysisTree/Cuts.hpp"
#include "CbmCuts.h"

#include "src/EntryConfig.hpp"
#include "src/Manager.hpp"
#include "src/Task.hpp"
#include "src/Utils.hpp"

const int nbins = 4000;
const float HugeValue = 1e9;

using namespace AnalysisTree;

struct Var
{
  std::string name_;
  float min_;
  float max_;
};
const std::vector<Var> momentum_ { {"px", -3, 3}, {"py", -3, 3}, {"pz", -2, 10},
                                  {"pT", 0, 3}, {"eta", -2, 8}, {"rapidity", 0, 5}, {"phi", -TMath::Pi(), TMath::Pi()} };
                                  
void LambdaCandidatesQA(QA::Task& task);

const std::string lambda_candidates_particles = "LambdaCandidates";
const std::string lambda_simulated_particles = "LambdaSimulated";

// SimpleCut nhitspos_cut({lambda_candidates_particles, "nhitspos"}, 8, 12);
// SimpleCut nhitsneg_cut({lambda_candidates_particles, "nhitsneg"}, 8, 12);
// SimpleCut sumnhits_cut({{lambda_candidates_particles, "nhitspos"}, {lambda_candidates_particles, "nhitsneg"}}, [](std::vector<double> par){ return par[0]+par[1] <= 14; });

// SimpleCut chi2primpos_cut({lambda_candidates_particles, "chi2primpos"}, 26, HugeValue);
// SimpleCut chi2primneg_cut({lambda_candidates_particles, "chi2primneg"}, 110, HugeValue);
// SimpleCut cosinepos_cut({lambda_candidates_particles, "cosinepos"}, 0.99825, HugeValue);
// SimpleCut ldl_cut({lambda_candidates_particles, "ldl"}, 4., HugeValue);

SimpleCut signal_cut({lambda_candidates_particles, "is_signal"}, 0, 2);


Cuts* selection_cuts = new Cuts("LambdaCandidatesCuts", {
//                                                           nhitspos_cut,
//                                                           nhitsneg_cut,
//                                                           sumnhits_cut,
//                                                           chi2primpos_cut,
//                                                           chi2primneg_cut,
//                                                           cosinepos_cut,
//                                                           ldl_cut,
                                                          signal_cut
                                                                          });

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./pfs_qa filelist" << std::endl;
    return -1;
  }

  const std::string filelist = argv[1];

//   QA::Manager man({filelist}, {"sTree"});
  QA::Manager man({filelist}, {"aTree"});
  man.SetOutFileName("pfsqa.root");
  
  man.AddBranchCut(selection_cuts);
  
  auto* task = new QA::Task;

  LambdaCandidatesQA(*task);
  
  man.AddTask(task);

  man.Init();
  man.Run(-1);
  man.Finish();

  return 0;
}

void LambdaCandidatesQA(QA::Task& task) {
  
//   Variable point_angle("PA", {{lambda_candidates_particles, "cosinepos"}}, [](std::vector<double> par){return TMath::ACos(par[0]);});
//   task.AddH1({"#alpha_{#Lambda, p}", point_angle, {nbins, 0, TMath::Pi()/3}});
// 
//   task.AddH1({"#chi^{2}_{prim, pos}", {lambda_candidates_particles, "chi2primpos"}, {nbins, 0, 400}});
//   task.AddH1({"#chi^{2}_{prim, neg}", {lambda_candidates_particles, "chi2primneg"}, {nbins, 0, 400}});
//   task.AddH1({"DCA, cm", {lambda_candidates_particles, "distance"}, {nbins, 0, 20}});
//   task.AddH1({"cos(#alpha_{#Lambda, p})", {lambda_candidates_particles, "cosinepos"}, {nbins, 0.5, 1}});
//   task.AddH1({"#chi^{2}_{geo}", {lambda_candidates_particles, "chi2geo"}, {nbins, 0, 100}});
//   task.AddH1({"L/#Delta L", {lambda_candidates_particles, "ldl"}, {nbins, 0, 100}});  
//   
//   task.AddH1({"N_{hits}, pos", {lambda_candidates_particles, "nhitspos"}, {10, 2.5, 12.5}});
//   task.AddH1({"N_{hits}, neg", {lambda_candidates_particles, "nhitsneg"}, {10, 2.5, 12.5}});
//   
//   Variable sumnhits("sumnhits", {{lambda_candidates_particles, "nhitspos"}, {lambda_candidates_particles, "nhitsneg"}}, [](std::vector<double> par){return par[0]+par[1];});
//   task.AddH1({"N_{hits}, sum", sumnhits, {19, 5.5, 24.5}});
  
  task.AddH1({"m_{#Lambda}, GeV", {lambda_candidates_particles, "mass"}, {nbins, 1, 2}});
  task.AddH1({"is signal", {lambda_candidates_particles, "is_signal"}, {5, -1.5, 3.5}});
/*  
  task.AddH2({"N_{hits}, pos", {lambda_candidates_particles, "nhitspos"}, {10, 2.5, 12.5}}, {"#chi^{2}_{prim, pos}", {lambda_candidates_particles, "chi2primpos"}, {nbins, 0, 400}});
  task.AddH2({"N_{hits}, neg", {lambda_candidates_particles, "nhitsneg"}, {10, 2.5, 12.5}}, {"#chi^{2}_{prim, neg}", {lambda_candidates_particles, "chi2primneg"}, {nbins, 0, 400}});
  task.AddH2({"N_{hits}, pos", sumnhits, {19, 5.5, 24.5}}, {"#chi^{2}_{geo}", {lambda_candidates_particles, "chi2geo"}, {nbins, 0, 100}});
  
  for(const auto& var : momentum_)
    task.AddH2({var.name_+", sim", {lambda_simulated_particles, var.name_}, {nbins, var.min_, var.max_}}, {var.name_+", rec", {lambda_candidates_particles, var.name_}, {nbins, var.min_, var.max_}});*/
}