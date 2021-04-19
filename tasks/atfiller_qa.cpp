#include "AnalysisTree/Cuts.hpp"
#include "CbmCuts.h"

#include "src/EntryConfig.hpp"
#include "src/Manager.hpp"
#include "src/Task.hpp"
#include "src/Utils.hpp"

using namespace AnalysisTree;

void ATFillerQA(QA::Task& task);

std::string eve_header = "AnaEventHeader";

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./atfiller_qa filelist1" << std::endl;
    return -1;
  }

  const std::string filelist1 = argv[1];
  
  QA::Manager man({filelist1}, {"cTree"});
  man.SetOutFileName("atfillerqa.root");
  
  auto* task = new QA::Task;

  ATFillerQA(*task);
  
  man.AddTask(task);

  man.Init();
  man.Run(-1);
  man.Finish();

  return 0;
}

void ATFillerQA(QA::Task& task) {
    
//   task.AddH1({"Multiplicity", {eve_header, "multiplicity"}, {1000, 0, 1000}});
  task.AddH1({"Centrality, %", {eve_header, "tracks_centrality"}, {20, 0, 100}});
}
