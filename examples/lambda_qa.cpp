#include "AnalysisTree/Cuts.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include "EntryConfig.hpp"
#include "Task.hpp"
#include "BasicQA.hpp"

using namespace AnalysisTree;

const std::string sim_particles = "SimParticles";
const std::string sim_event_header = "SimEventHeader";

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./lambda_qa filelist" << std::endl;
    return -1;
  }
  int n_events = -1;
  
  if (argc == 3) {
    n_events = atoi(argv[2]);
  }

  const std::string filelist = argv[1];

  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("lambda_qa.root");

  auto centr_cut = RangeCut("SimEventHeader.b", 0, 3);
  auto lambda_cut = EqualsCut("SimParticles.pid", 3122);
  auto prim_cut = EqualsCut("SimParticles.mother_id", -1);

  Cuts* centr = new Cuts("SemiCentral", {centr_cut});
  task->SetEventCuts(centr);
  
  Cuts* mc_lambda_primiry = new Cuts("McLambda", {lambda_cut, prim_cut});
  task->AddH1({"b (fm)", {sim_event_header, "b"}, {100, 0, 20}},mc_lambda_primiry);
  task->AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {15, 0., 3.}}, {"p_{T}, GeV/c", {sim_particles, "pT"}, {15, 0., 3.}}, mc_lambda_primiry);
  AddParticleQA(task, sim_particles, mc_lambda_primiry);
  
  
  
  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->Run(n_events);
  man->Finish();

  return 0;
}


