#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 2100;

using namespace AnalysisTree;

void pfs_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("tof_purities_v0_qa.root");

  struct Particle
  {
    std::string name_;
    int pdg_;
    std::string daughter_neg_name_;
    std::string daughter_pos_name_;  
  };
  
  std::vector<Particle> particles
  {
    {"Lambda", 3122, "pi", "p"},
    {"Kshort", 310, "pi", "pi"}
  };
  
  for(auto& pa : particles)
  {
    SimpleCut mother_pdg_cut = EqualsCut("RecParticles.mother_pdg", pa.pdg_);
    
    SimpleCut q_neg_cut = EqualsCut("RecParticles.q", -1);
    Cuts* selection_neg_daughter = new Cuts((pa.name_ + "_" + pa.daughter_neg_name_ + "_minus_sgnl").c_str(), {mother_pdg_cut, q_neg_cut});
    task->AddH1({"purity", Variable::FromString(("RecParticles.prob_" + pa.daughter_neg_name_).c_str()), {nbins, -1.05, 1.05}}, selection_neg_daughter);
    
    SimpleCut q_pos_cut = EqualsCut("RecParticles.q", +1);
    Cuts* selection_pos_daughter = new Cuts((pa.name_ + "_" +  pa.daughter_pos_name_ + "_plus_sgnl").c_str(), {mother_pdg_cut, q_pos_cut});
    task->AddH1({"purity", Variable::FromString(("RecParticles.prob_" + pa.daughter_pos_name_).c_str()), {nbins, -1.05, 1.05}}, selection_pos_daughter);
  }
  
  man->AddTask(task);

  man->Init({filelist}, {"aTree"});
  man->Run(-1);
  man->Finish();
  
}

int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./pfs_qa filelist" << std::endl;
    return -1;
  }  
  
  pfs_qa(argv[1]);
  return 0;
}