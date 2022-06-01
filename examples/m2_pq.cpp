#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 1000;

using namespace AnalysisTree;

void pfs_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("m2_pq.root");
  
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
    {"bg", 1, true},
    {"none", 2, true}
  };
  
  for(auto& pa : particles)
  {
//     double pdg = pa.pdg_;
    
    SimpleCut pdg_simplecut = EqualsCut("RecParticles.pid", pa.pdg_);
//     SimpleCut pdg_simplecut = RangeCut("RecParticles.pid", pa.pdg_-0.1, pa.pdg_+0.1);
//     SimpleCut pdg_simplecut = SimpleCut({"RecParticles.pid"}, [pdg]( std::vector<double>& var ) { return std::fabs(var.at(0)-pdg)<0.1;});
    Cuts* pdg_cut = new Cuts((pa.name_ + "_plus").c_str(), {pdg_simplecut});
    task->AddH2({"p/q, GeV/c", Variable::FromString("TofHits.qp_tof"), {nbins, -10, 10}}, {"m^{2}/q^{2}, (GeV/c^{2})^{2}", Variable::FromString("TofHits.mass2"), {nbins, -2, 5}}, pdg_cut);
    
    if(!pa.is_antiparticle_) continue;
    
    SimpleCut antipdg_simplecut = EqualsCut("RecParticles.pid", -pa.pdg_);
    Cuts* antipdg_cut = new Cuts((pa.name_ + "_minus").c_str(), {antipdg_simplecut});
    task->AddH2({"p/q, GeV/c", Variable::FromString("TofHits.qp_tof"), {nbins, -10, 10}}, {"m^{2}/q^{2}, (GeV/c^{2})^{2}", Variable::FromString("TofHits.mass2"), {nbins, -2, 5}}, antipdg_cut);
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