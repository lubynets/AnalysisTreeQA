#include "AnalysisTree/Cuts.hpp"
#include "CbmCuts.h"

#include "src/EntryConfig.hpp"
#include "src/Manager.hpp"
#include "src/Task.hpp"
#include "src/Utils.hpp"

const int nbins = 600;
const float HugeValue = 1e9;

using namespace AnalysisTree;

void mass3D(QA::Task& task);
std::string StringBinNumber(int number);

const std::string reco_tracks = "RecParticlesMcPid";
const std::string event_header = "AnaEventHeader";
const float y_beam = 1.62179;

// SimpleCut signal_cut({reco_tracks, "is_signal"}, 0, 1);

// Cuts* selection_cuts = new Cuts("LambdaCandidatesCuts", {signal_cut});

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./mass3D filelist" << std::endl;
    return -1;
  }

  const std::string filelist = argv[1];

  QA::Manager man({filelist}, {"cTree"});
  man.SetOutFileName("out.mass3D.root");
  
//   man.AddBranchCut(selection_cuts);
  
  auto* task = new QA::Task;

  mass3D(*task);
  
  man.AddTask(task);

  man.Init();
  man.Run(-1);
  man.Finish();

  return 0;
}

void mass3D(QA::Task& task) {
  
  const int y_nbins = 5;
  const float y_low = y_beam-0.5;
  const float y_up = y_beam+1.0;
  const float y_bin_width = (y_up-y_low)/y_nbins;

  const int pT_nbins = 5;
  const float pT_low = 0.3;
  const float pT_up = 1.3;
  const float pT_bin_width = (pT_up-pT_low)/pT_nbins;
  
  const std::vector<float> C_binranges {0, 20, 40, 100};
  const int C_nbins = C_binranges.size()-1;
  
  std::vector<Cuts*> cut_3D;
  cut_3D.resize(y_nbins*pT_nbins*C_nbins);
  
  for(int i_C_bin=1; i_C_bin<=C_nbins; i_C_bin++)
    for(int i_y_bin=1; i_y_bin<=y_nbins; i_y_bin++)
      for(int i_pT_bin=1; i_pT_bin<=pT_nbins; i_pT_bin++)
      {
        std::string binname = "C" + StringBinNumber(i_C_bin) + "_y" + StringBinNumber(i_y_bin) + "_pT" + StringBinNumber(i_pT_bin);
        int binnumber = y_nbins*pT_nbins*(i_C_bin-1) + pT_nbins*(i_y_bin-1) + (i_pT_bin-1);
              
        SimpleCut y_cut({reco_tracks, "rapidity"}, y_low + (i_y_bin-1)*y_bin_width, y_low + i_y_bin*y_bin_width);
        SimpleCut pT_cut({reco_tracks, "pT"}, pT_low + (i_pT_bin-1)*pT_bin_width, pT_low + i_pT_bin*pT_bin_width);
        SimpleCut C_cut({event_header, "tracks_centrality"}, C_binranges.at(i_C_bin-1), C_binranges.at(i_C_bin));
        
        cut_3D.at(binnumber) = new Cuts(binname, {C_cut, y_cut, pT_cut});
        
        task.AddH1({"m_{#Lambda}, GeV/c^{2}", {reco_tracks, "mass"}, {nbins, 1.0, 1.232}}, cut_3D.at(binnumber));
      }
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}
