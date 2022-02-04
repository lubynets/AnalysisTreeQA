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
  task->SetOutputFileName("tof_purities_qa.root");
  
  SimpleCut pdgs_cut = SimpleCut({"RecParticles.mc_pdg"}, []( std::vector<double>& var ) { return std::fabs(var.at(0)-211)<0.1 || 
                                                                                                  std::fabs(var.at(0)+211)<0.1 ||
                                                                                                  std::fabs(var.at(0)-321)<0.1 ||
                                                                                                  std::fabs(var.at(0)+321)<0.1 ||
                                                                                                  std::fabs(var.at(0)-2212)<0.1 ||
                                                                                                  std::fabs(var.at(0)-1000010020)<0.1;});
  
  SimpleCut qs_cut = SimpleCut({"RecParticles.q"}, []( std::vector<double>& var ) { return std::fabs(var.at(0)-1)<0.1 || 
                                                                                           std::fabs(var.at(0)+1)<0.1;});
  
  struct Particle
  {
    std::string name_;
    int pdg_;
    bool is_antiparticle_;
  };
  
  std::vector<Particle> particles
  {
    {"pi", 211, true},
    {"K", 321, true},
    {"p", 2212, false},
    {"d", 1000010020, false},
//     {"he3", 1000020030, false},
  };
  
  for(auto& pa : particles)
  {
    std::vector<SimpleCut> current_cuts;
    current_cuts.push_back(pdgs_cut);
    current_cuts.push_back(qs_cut);
    
    double pdg = pa.pdg_;
    
    SimpleCut sgnl_cut_plus = SimpleCut({"RecParticles.mc_pdg", "RecParticles.q"}, [pdg]( std::vector<double>& var ) { return std::fabs(var.at(0)-pdg)<0.1 && std::fabs(var.at(1)-1)<0.1;});
    Cuts* selection_sgnl_plus = new Cuts((pa.name_ + "_plus_sgnl").c_str(), {pdgs_cut, qs_cut, sgnl_cut_plus});
    task->AddH1({"purity", Variable::FromString(("RecParticles.prob_" + pa.name_).c_str()), {nbins, -1.05, 1.05}}, selection_sgnl_plus);
    
    SimpleCut bckgr_cut_plus = SimpleCut({"RecParticles.mc_pdg", "RecParticles.q"}, [pdg]( std::vector<double>& var ) { return std::fabs(var.at(0)-pdg)>0.1 && std::fabs(var.at(1)-1)<0.1;});
    Cuts* selection_bckgr_plus = new Cuts((pa.name_ + "_plus_bckgr").c_str(), {pdgs_cut, qs_cut, bckgr_cut_plus});
    task->AddH1({"purity", Variable::FromString(("RecParticles.prob_" + pa.name_).c_str()), {nbins, -1.05, 1.05}}, selection_bckgr_plus);
    
    if(!pa.is_antiparticle_) continue;
    
    SimpleCut sgnl_cut_minus = SimpleCut({"RecParticles.mc_pdg", "RecParticles.q"}, [pdg]( std::vector<double>& var ) { return std::fabs(var.at(0)+pdg)<0.1 && std::fabs(var.at(1)+1)<0.1;});
    Cuts* selection_sgnl_minus = new Cuts((pa.name_ + "_minus_sgnl").c_str(), {pdgs_cut, qs_cut, sgnl_cut_minus});
    task->AddH1({"purity", Variable::FromString(("RecParticles.prob_" + pa.name_).c_str()), {nbins, -1.05, 1.05}}, selection_sgnl_minus);
    
    SimpleCut bckgr_cut_minus = SimpleCut({"RecParticles.mc_pdg", "RecParticles.q"}, [pdg]( std::vector<double>& var ) { return std::fabs(var.at(0)+pdg)>0.1 && std::fabs(var.at(1)+1)<0.1;});
    Cuts* selection_bckgr_minus = new Cuts((pa.name_ + "_minus_bckgr").c_str(), {pdgs_cut, qs_cut, bckgr_cut_minus});
    task->AddH1({"purity", Variable::FromString(("RecParticles.prob_" + pa.name_).c_str()), {nbins, -1.05, 1.05}}, selection_bckgr_minus);
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