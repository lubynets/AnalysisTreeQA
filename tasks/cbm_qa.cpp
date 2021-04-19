#include "AnalysisTree/Cuts.hpp"
#include "CbmCuts.h"

#include "src/EntryConfig.hpp"
#include "src/Manager.hpp"
#include "src/Task.hpp"
#include "src/Utils.hpp"

typedef AnalysisTree::QA::EntryConfig::PlotType PlotType;
using AnalysisTree::QA::gNbins;

//namespace AnalysisTree{
//namespace QA{

using namespace AnalysisTree;

void KFPFTracksQA(QA::Task& task);
void VertexTracksQA(QA::Task& task);
void TofHitsQA(QA::Task& task);
void SimParticlesQA(QA::Task& task);
void SimEventHeaderQA(QA::Task& task);
void RecEventHeaderQA(QA::Task& task);
void EfficiencyMaps(QA::Task& task);

const std::string sim_event_header = "SimEventHeader";
const std::string rec_event_header = "RecEventHeader";
const std::string sim_particles = "SimParticles";
const std::string rec_tracks = "VtxTracks";
const std::string tof_hits = "TofHits";

std::vector<size_t> pdf_plots{};

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./cbm_qa filelist" << std::endl;
    return -1;
  }

  const std::string filelist = argv[1];

  QA::Manager man({filelist}, {"aTree"});
  man.SetOutFileName("cbmqa.root");

  man.AddBranchCut(GetCbmTrackCuts(rec_tracks));
  man.AddBranchCut(GetCbmMcTracksCuts(sim_particles));
  man.AddBranchCut(GetCbmTofHitsCuts(tof_hits));
  man.SetEventCuts(GetCbmEventCuts(rec_event_header));

  auto* task = new QA::Task;

  KFPFTracksQA(*task);
  VertexTracksQA(*task);
  TofHitsQA(*task);
  SimParticlesQA(*task);
  SimEventHeaderQA(*task);
  RecEventHeaderQA(*task);
  EfficiencyMaps(*task);

  AddParticlesFlowQA(task, sim_particles, {sim_event_header, "psi_RP"}, {2212, 211, -211});

  man.AddTask(task);

  man.Init();
  man.Run(-1);
  man.Finish();

  return 0;
}

void VertexTracksQA(QA::Task& task) {
  AddTrackQA(&task, rec_tracks);
  if (!sim_particles.empty()) {
    AddTracksMatchQA(&task, rec_tracks, sim_particles);
  }

  Variable chi2_over_ndf("chi2_ndf", {{rec_tracks, "chi2"}, {rec_tracks, "ndf"}}, [](std::vector<double>& var) { return var.at(0) / var.at(1); });

  task.AddH1({"DCA_{x}, cm", {rec_tracks, "dcax"}, {gNbins, -1, 1}});
  task.AddH1({"DCA_{y}, cm", {rec_tracks, "dcay"}, {gNbins, -1, 1}});
  task.AddH1({"DCA_{z}, cm", {rec_tracks, "dcaz"}, {gNbins, -1, 1}});
  task.AddH1({"NDF", {rec_tracks, "ndf"}, {30, 0, 30}});
  task.AddH1({"Number of hits (STS+MVD)", {rec_tracks, "nhits"}, {15, 0, 15}});
  task.AddH1({"Number of hits (MVD)", {rec_tracks, "nhits_mvd"}, {5, 0, 5}});
  task.AddH1({"#chi^{2}_{vertex}", {rec_tracks, "vtx_chi2"}, {500, 0, 100}});
  task.AddH1({"#chi^{2}/NDF", chi2_over_ndf, {gNbins, 0, 10}});
  task.AddH2({"DCA_{x}, cm", {rec_tracks, "dcax"}, {gNbins, -1, 1}}, {"DCA_{y}, cm", {rec_tracks, "dcay"}, {gNbins, -1, 1}});
}

void TofHitsQA(QA::Task& task) {
  task.AddH1({"TOF hit x-position (cm)", {tof_hits, "x"}, {gNbins, -600, 600}});
  task.AddH1({"TOF hit y-position (cm)", {tof_hits, "y"}, {gNbins, -400, 400}});
  task.AddH1({"TOF hit z-position (cm)", {tof_hits, "z"}, {gNbins, 600, 800}});

  task.AddH2({"TOF hit x-position (cm)", {tof_hits, "x"}, {gNbins, -600, 600}}, {"TOF hit y-position (cm)", {tof_hits, "y"}, {gNbins, -600, 600}});

  Variable qp_sts("qp_reco", {{rec_tracks, "q"}, {rec_tracks, "p"}}, [](std::vector<double>& qp) { return qp.at(0) * qp.at(1); });
  task.AddH2({"sign(q)*p, GeV/c", qp_sts, {gNbins, -10, 10}}, {"m^{2}, GeV^{2}/c^{2}", {tof_hits, "mass2"}, {gNbins, -1, 2}});
  task.AddH2({"sign(q)*p STS, GeV/c", qp_sts, {gNbins, -10, 10}}, {"sign(q)*p TOF, GeV/c", {tof_hits, "qp_tof"}, {gNbins, -10, 10}});
}

void SimParticlesQA(QA::Task& task) {
  AddParticleQA(&task, sim_particles);

  task.AddH1({"PDG code", {sim_particles, "pid"}, {8001, -3999.5, 4000.5}});
}

void SimEventHeaderQA(QA::Task& task) {
  task.AddH1({"x_{vertex}^{MC} (cm)", {sim_event_header, "vtx_x"}, {gNbins, -1, 1}});
  task.AddH1({"y_{vertex}^{MC} (cm)", {sim_event_header, "vtx_y"}, {gNbins, -1, 1}});
  task.AddH1({"z_{vertex}^{MC} (cm)", {sim_event_header, "vtx_z"}, {gNbins, -1, 1}});
  task.AddH1({"b (fm)", {sim_event_header, "b"}, {gNbins, 0, 20}});
  task.AddH1({"#Psi_{RP}", {sim_event_header, "psi_RP"}, {gNbins, 0, 6.5}});

  task.AddH2({"x_{vertex}^{MC} (cm)", {sim_event_header, "vtx_x"}, {gNbins, -1, 1}}, {"y_{vertex}^{MC} (cm)", {sim_event_header, "vtx_y"}, {gNbins, -1, 1}});
}

void RecEventHeaderQA(QA::Task& task) {
  task.AddH1({"x_{vertex} (cm)", {rec_event_header, "vtx_x"}, {gNbins, -1, 1}});
  task.AddH1({"y_{vertex} (cm)", {rec_event_header, "vtx_y"}, {gNbins, -1, 1}});
  task.AddH1({"z_{vertex} (cm)", {rec_event_header, "vtx_z"}, {gNbins, -1, 1}});
  task.AddH1({"#chi^{2}_{vertex fit}", {rec_event_header, "vtx_chi2"}, {gNbins, 0, 5}});

  task.AddH1({"E_{PSD} (GeV)", {rec_event_header, "Epsd"}, {gNbins, 0, 60}});
  task.AddH1({"M_{tracks}", {rec_event_header, "M"}, {800, 0, 800}});
  task.AddH1({"Event ID", {rec_event_header, "evt_id"}, {gNbins, 0, 2000}});

  task.AddH2({"x_{vertex} (cm)", {rec_event_header, "vtx_x"}, {gNbins, -1, 1}}, {"y_{vertex} (cm)", {rec_event_header, "vtx_y"}, {gNbins, -1, 1}});

  task.AddH2({"x_{vertex} (cm)", {rec_event_header, "vtx_x"}, {gNbins, -1, 1}}, {"x_{vertex} (cm)", {sim_event_header, "vtx_x"}, {gNbins, -1, 1}});
  task.AddH2({"y_{vertex} (cm)", {rec_event_header, "vtx_y"}, {gNbins, -1, 1}}, {"y_{vertex} (cm)", {sim_event_header, "vtx_y"}, {gNbins, -1, 1}});
  task.AddH2({"z_{vertex} (cm)", {rec_event_header, "vtx_z"}, {gNbins, -1, 1}}, {"z_{vertex} (cm)", {sim_event_header, "vtx_z"}, {gNbins, -1, 1}});
}

void EfficiencyMaps(QA::Task& task) {

  const float y_beam = 1.62179f;// TODO from DataHeader
  const float p_mass = AnalysisTree::GetMassByPdgId(2212);
  const float pi_mass = AnalysisTree::GetMassByPdgId(211);

  Variable proton_y("y-y_{beam}", {{rec_tracks, "p"}, {rec_tracks, "pz"}},
                    [y_beam, p_mass](std::vector<double>& var) {
                        const float e = sqrt(p_mass*p_mass + var[0]*var[0]);
                        return 0.5*log( (e+var[1]) / (e-var[1]) ); });

  Variable pion_y("y-y_{beam}", {{rec_tracks, "p"}, {rec_tracks, "pz"}},
                  [y_beam, pi_mass](std::vector<double>& var) {
                      const float e = sqrt(pi_mass*pi_mass + var[0]*var[0]);
                      return 0.5*log( (e+var[1]) / (e-var[1]) ); });

  Cuts* mc_protons = new Cuts("McProtons", {{{sim_particles, "pid"}, 2212}});
  Cuts* mc_pions_neg = new Cuts("McPionsNeg", {{{sim_particles, "pid"}, -211}});
  Cuts* mc_pions_pos = new Cuts("McPionsPos", {{{sim_particles, "pid"}, 211}});

  task.AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {gNbins, -1, 5}}, {"p_{T}, GeV/c", {sim_particles, "pT"}, {gNbins, 0, 2}}, mc_protons);
  task.AddH2({"#it{y}_{Lab}", proton_y, {gNbins, -1, 5}}, {"p_{T}, GeV/c", {rec_tracks, "pT"}, {gNbins, 0, 2}}, mc_protons);

  task.AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {gNbins, -1, 5}}, {"p_{T}, GeV/c", {sim_particles, "pT"}, {gNbins, 0, 2}}, mc_pions_neg);
  task.AddH2({"#it{y}_{Lab}", pion_y, {gNbins, -1, 5}}, {"p_{T}, GeV/c", {rec_tracks, "pT"}, {gNbins, 0, 2}}, mc_pions_neg);

  task.AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {gNbins, -1, 5}}, {"p_{T}, GeV/c", {sim_particles, "pT"}, {gNbins, 0, 2}}, mc_pions_pos);
  task.AddH2({"#it{y}_{Lab}", pion_y, {gNbins, -1, 5}}, {"p_{T}, GeV/c", {rec_tracks, "pT"}, {gNbins, 0, 2}}, mc_pions_pos);
}

void KFPFTracksQA(QA::Task& task) {
  const std::vector<std::string> field_par_names_{"cx0", "cx1", "cx2", "cy1", "cy2", "cz0", "cz1", "cz2"};
  const std::vector<std::string> cov_names_{"cov1", "cov2", "cov3", "cov4", "cov5", "cov6", "cov7", "cov8", "cov9",
                                            "cov10", "cov11", "cov12", "cov13", "cov14", "cov15"};
  // KFPF tracks
  task.AddH1({"x, cm", {rec_tracks, "x"}, {gNbins, -10, 10}});
  task.AddH1({"y, cm", {rec_tracks, "y"}, {gNbins, -10, 10}});
  task.AddH1({"z, cm", {rec_tracks, "z"}, {gNbins, 0, 40}});
  task.AddH1({"t_{x}", {rec_tracks, "tx"}, {gNbins, -2, 2}});
  task.AddH1({"t_{y}", {rec_tracks, "ty"}, {gNbins, -2, 2}});
  task.AddH1({"q*p, GeV/c", {rec_tracks, "qp"}, {gNbins, -10, 10}});

  for (const auto& field_par_name : field_par_names_)
    task.AddH1({field_par_name, {rec_tracks, field_par_name}, {gNbins, -1, 1}});

  for (const auto& cov_name : cov_names_)
    task.AddH1({cov_name, {rec_tracks, cov_name}, {1000, -.1, .1}});

  task.AddH1({"cy0", {rec_tracks, "cy0"}, {gNbins, -12, -8}});
  task.AddH1({"z0", {rec_tracks, "z0"}, {gNbins, 0, 40}});
}

//}
//}

//void SavePdf(QA::Task* task){
//
//  Helper::Canvas t("");
//  t.SetCanvasSize(1200, 900);
//  t.SetLabelPosition(0.5, 0.87);
//  t.AddLabel("CBM simulation, GEANT3");
//  t.AddLabel("UrQMD, EoS=0");
//  t.AddLabel("Au+Au @ 12#it{A} GeV/#it{c}");
//
//  for(auto plot_id : pdf_plots){
//    auto plot = task->Entries().at(plot_id);
//
//    switch(plot.GetType()) {
//      case PlotType::kHisto1D : t.Add1DPlot(plot.GetPlot() )
//
//    }
//  }
//
//}
