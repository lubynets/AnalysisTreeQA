#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "BasicQA.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void example(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("cbm_qa_lite.root");
  
  struct Particle {
    std::string name_;
    SimpleCut cut_;
  };
  
  std::vector<Particle> particles {
    {"kaonplus", EqualsCut("RecParticles.pid", 321)},
    {"kaonminus", EqualsCut("RecParticles.pid", -321)},
    {"pionplus", EqualsCut("RecParticles.pid", 211)},
    {"pionminus", EqualsCut("RecParticles.pid", -211)},
    {"proton", EqualsCut("RecParticles.pid", 2212)}
  };
  
  for(auto& pa : particles) {
    Cuts* recpid_cut = new Cuts(pa.name_.c_str(), {pa.cut_});
    AddParticleQA(task, "RecParticles", recpid_cut);
  }
   
  man->AddTask(task);

  man->Init({filelist}, {"aTree"});
  man->Run(-1);
  man->Finish();
}


int main(int argc, char* argv[]){
    if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./cbm_qa_lite filename" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filename = argv[1];
  example(filename);

  return 0;
}  