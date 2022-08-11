#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void example(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("v1_qa.root");
  
  SimpleCut mc_primary =  EqualsCut("SimParticles.mother_id", -1);
  
  struct Particle {
    std::string name_;
    SimpleCut cut_;
  };
  
  std::vector<Particle> particles {
    {"kaonplus", EqualsCut("SimParticles.pid", 321)},
    {"kaonminus", EqualsCut("SimParticles.pid", -321)},
    {"kshort", EqualsCut("SimParticles.pid", 310)}
  };
  
  struct Pt_Range {
    std::string name_;
    SimpleCut cut_;
  };
  
  std::vector<Pt_Range> pt_ranges {
    {"0_0__0_4", RangeCut("SimParticles.pT", 0, 0.4)},
    {"0_4__0_8", RangeCut("SimParticles.pT", 0.4, 0.8)},
    {"0_8__1_2", RangeCut("SimParticles.pT", 0.8, 1.2)},
    {"0_3__1_35", RangeCut("SimParticles.pT", 0.3, 1.35)}
  };
    
  const Field psi_RP = Field("SimEventHeader", "psi_RP");
  const Field mc_phi = Field("SimParticles", "phi");
  Variable v1("v1", {mc_phi, psi_RP}, [](std::vector<double> phi) { return TMath::Cos(phi[0] - phi[1]); });
//   Variable v2("v2", {mc_phi, psi_RP}, [](std::vector<double> phi) { return TMath::Cos(2*(phi[0] - phi[1])); });
  
  std::vector<Cuts*> current_cuts;
  current_cuts.resize(particles.size()*pt_ranges.size());
  int binnumber{0};
    
  for(auto& pa : particles) {
    for(auto& pr : pt_ranges) {
      std::string cuts_name = pa.name_ + "_" + pr.name_;
      current_cuts.at(binnumber) = new Cuts(cuts_name.c_str(), {mc_primary, pa.cut_, pr.cut_});
      task->AddProfile({"#it{y}", Variable::FromString("SimParticles.rapidity"), {11, 0.52179, 2.72179}}, {"v_{1}", v1, {}}, current_cuts.at(binnumber));
      binnumber++;
    }
  }
  
  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->Run(-1);
  man->Finish();
}

int main(int argc, char* argv[]){
    if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./v1_qa filename" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filename = argv[1];
  example(filename);

  return 0;
}  