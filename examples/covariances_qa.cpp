#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

const int nbins = 1000;

using namespace AnalysisTree;

void covariances_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("cbm_qa.root");
  
  SimpleCut passcuts = EqualsCut("VtxTracks.pass_cuts", true);
  Cuts* cuts = new Cuts("passcuts", {passcuts});
  
  std::vector<std::string> covariances {"cov1", "cov2", "cov3", "cov4", "cov5", "cov6", "cov7", "cov8", "cov9", "cov10",
                                        "cov11", "cov12", "cov13", "cov14", "cov15"};
  
  Variable chi2_over_ndf("chi2ndf", {{"VtxTracks", "chi2"}, {"VtxTracks", "ndf"}}, []( std::vector<double>& var ) { return var.at(0)/var.at(1); });
  task->AddH1({"#chi^{2}/NDF", chi2_over_ndf, {nbins, 0, 20}}, cuts);                               
                                        
  for(auto& cov : covariances)  
  {
    task->AddH1({cov.c_str(), Variable::FromString(("VtxTracks." + cov).c_str()), {nbins, -0.0001, 0.0001}}, cuts);
    task->AddH2({cov.c_str(), Variable::FromString(("VtxTracks." + cov).c_str()), {nbins, -0.0001, 0.0001}}, {"#chi^{2}/NDF", chi2_over_ndf, {nbins, 0, 20}}, cuts);
    task->AddProfile({cov.c_str(), Variable::FromString(("VtxTracks." + cov).c_str()), {nbins, -0.0001, 0.0001}}, {"#chi^{2}/NDF", chi2_over_ndf, {nbins, 0, 20}}, cuts);
    task->AddProfile({"#chi^{2}/NDF", chi2_over_ndf, {nbins, 0, 20}}, {cov.c_str(), Variable::FromString(("VtxTracks." + cov).c_str()), {nbins, -0.0001, 0.0001}}, cuts);
  }
    

  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->Run(-1);
  man->Finish();  
}

int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./covariances_qa filelist" << std::endl;
    return -1;
  }  
  
  covariances_qa(argv[1]);
  return 0;
}