#include "AnalysisTree/Cuts.hpp"
#include "CbmCuts.h"

#include "src/EntryConfig.hpp"
#include "src/Manager.hpp"
#include "src/Task.hpp"
#include "src/Utils.hpp"

const int nbins = 500;
const float HugeValue = 1e9;

using namespace AnalysisTree;

void Qa3D(QA::Task& task);

const std::string lambda_candidates_particles = "LambdaCandidates";
const std::string lambda_simulated_particles = "LambdaSimulated";
const float y_beam = 1.62179;

SimpleCut signal_cut({lambda_candidates_particles, "is_signal"}, 0, 1);

Cuts* selection_cuts = new Cuts("LambdaCandidatesCuts", {signal_cut});

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./qa3D filelist" << std::endl;
    return -1;
  }

  const std::string filelist = argv[1];

  QA::Manager man({filelist}, {"sTree"});
  man.SetOutFileName("qa3d.root");
  
  man.AddBranchCut(selection_cuts);
  
  auto* task = new QA::Task;

  Qa3D(*task);
  
  man.AddTask(task);

  man.Init();
  man.Run(-1);
  man.Finish();

  return 0;
}

void Qa3D(QA::Task& task) {
  
  const int y_nbins = 20;
  const float y_low = y_beam-0.8;
  const float y_up = y_beam+1.2;
  const float y_bin_width = (y_up-y_low)/y_nbins;

  const int pT_nbins = 18;
  const float pT_low = 0.;
  const float pT_up = 1.8;
  const float pT_bin_width = (pT_up-pT_low)/pT_nbins;
  
  std::vector<Cuts*> y_pT_cut;
  y_pT_cut.resize(y_nbins*pT_nbins);
  
  Variable diff_px("diff_px", {{lambda_candidates_particles, "px"}, {lambda_simulated_particles, "px"}}, [](std::vector<double> par){return par[0]-par[1];});
  Variable diff_py("diff_py", {{lambda_candidates_particles, "py"}, {lambda_simulated_particles, "py"}}, [](std::vector<double> par){return par[0]-par[1];});
  Variable diff_pz("diff_pz", {{lambda_candidates_particles, "pz"}, {lambda_simulated_particles, "pz"}}, [](std::vector<double> par){return par[0]-par[1];});  
        
  Variable pull_px("pull_px", {{lambda_candidates_particles, "px"}, {lambda_simulated_particles, "px"}, {lambda_candidates_particles, "pxerr"}}, [](std::vector<double> par){return (par[0]-par[1])/par[2];});
  Variable pull_py("pull_py", {{lambda_candidates_particles, "py"}, {lambda_simulated_particles, "py"}, {lambda_candidates_particles, "pyerr"}}, [](std::vector<double> par){return (par[0]-par[1])/par[2];});
  Variable pull_pz("pull_pz", {{lambda_candidates_particles, "pz"}, {lambda_simulated_particles, "pz"}, {lambda_candidates_particles, "pzerr"}}, [](std::vector<double> par){return (par[0]-par[1])/par[2];});  
        
  for(int i_y_bin=1; i_y_bin<=y_nbins; i_y_bin++)
    for(int i_pT_bin=1; i_pT_bin<=pT_nbins; i_pT_bin++)
    {
      std::string binname;
      if(i_y_bin<10 && i_pT_bin<10)
        binname = "y0" + std::to_string(i_y_bin) + "_pT0" + std::to_string(i_pT_bin);
      else if(i_y_bin<10)
        binname = "y0" + std::to_string(i_y_bin) + "_pT" + std::to_string(i_pT_bin);
      else if(i_pT_bin<10)
        binname = "y" + std::to_string(i_y_bin) + "_pT0" + std::to_string(i_pT_bin);
      else      
        binname = "y" + std::to_string(i_y_bin) + "_pT" + std::to_string(i_pT_bin);
      
      int binnumber = (i_y_bin-1)*pT_nbins + i_pT_bin-1;
            
      SimpleCut y_cut({lambda_candidates_particles, "rapidity"}, y_low + (i_y_bin-1)*y_bin_width, y_low + i_y_bin*y_bin_width);
      SimpleCut pT_cut({lambda_candidates_particles, "pT"}, pT_low + (i_pT_bin-1)*pT_bin_width, pT_low + i_pT_bin*pT_bin_width);
      
      y_pT_cut.at(binnumber) = new Cuts(binname, {y_cut, pT_cut});
      
      task.AddH1({"m_{#Lambda}, GeV/c^{2}", {lambda_candidates_particles, "mass"}, {nbins, 1.1, 1.132}}, y_pT_cut.at(binnumber));
      
      task.AddH1({"p_{x}^{reco}-p_{x}^{sim}, GeV/c", diff_px, {nbins, -0.2, 0.2}}, y_pT_cut.at(binnumber));
      task.AddH1({"p_{y}^{reco}-p_{y}^{sim}, GeV/c", diff_py, {nbins, -0.2, 0.2}}, y_pT_cut.at(binnumber));
      task.AddH1({"p_{z}^{reco}-p_{z}^{sim}, GeV/c", diff_pz, {nbins, -0.5, 0.5}}, y_pT_cut.at(binnumber));
      
      task.AddH2({"p_{x}^{reco}-p_{x}^{sim}, GeV/c", diff_px, {nbins, -0.2, 0.2}}, {"p_{y}^{reco}-p_{y}^{sim}, GeV/c", diff_py, {nbins, -0.2, 0.2}}, y_pT_cut.at(binnumber));
      task.AddH2({"p_{x}^{reco}-p_{x}^{sim}, GeV/c", diff_px, {nbins, -0.2, 0.2}}, {"p_{z}^{reco}-p_{z}^{sim}, GeV/c", diff_pz, {nbins, -0.5, 0.5}}, y_pT_cut.at(binnumber));
      task.AddH2({"p_{y}^{reco}-p_{y}^{sim}, GeV/c", diff_py, {nbins, -0.2, 0.2}}, {"p_{z}^{reco}-p_{z}^{sim}, GeV/c", diff_pz, {nbins, -0.5, 0.5}}, y_pT_cut.at(binnumber));
      
      task.AddH1({"(p_{x}^{reco}-p_{x}^{sim})/#sigma_{p_{x}}", pull_px, {nbins, -10, 10}}, y_pT_cut.at(binnumber));
      task.AddH1({"(p_{y}^{reco}-p_{y}^{sim})/#sigma_{p_{y}}", pull_py, {nbins, -10, 10}}, y_pT_cut.at(binnumber));
      task.AddH1({"(p_{z}^{reco}-p_{z}^{sim})/#sigma_{p_{z}}", pull_pz, {nbins, -10, 10}}, y_pT_cut.at(binnumber));
       
      task.AddH2({"(p_{x}^{reco}-p_{x}^{sim})/#sigma_{p_{x}}", pull_px, {nbins, -10, 10}}, {"(p_{y}^{reco}-p_{y}^{sim})/#sigma_{p_{y}}", pull_py, {nbins, -10, 10}}, y_pT_cut.at(binnumber));
      task.AddH2({"(p_{x}^{reco}-p_{x}^{sim})/#sigma_{p_{x}}", pull_px, {nbins, -10, 10}}, {"(p_{z}^{reco}-p_{z}^{sim})/#sigma_{p_{z}}", pull_pz, {nbins, -10, 10}}, y_pT_cut.at(binnumber));
      task.AddH2({"(p_{y}^{reco}-p_{y}^{sim})/#sigma_{p_{y}}", pull_py, {nbins, -10, 10}}, {"(p_{z}^{reco}-p_{z}^{sim})/#sigma_{p_{z}}", pull_pz, {nbins, -10, 10}}, y_pT_cut.at(binnumber));
      
    }
}