#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 1000;
const float HugeValue = 1e9;
const float PI = TMath::Pi();

using namespace AnalysisTree;

void pfs_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("pfs_qa.root");
  
  SimpleCut allcut = RangeCut("Candidates.generation", -0.1, 12.1);
  SimpleCut sgnlcut = RangeCut("Candidates.generation",  0.9, 12.1);
  SimpleCut bckgrcut = RangeCut("Candidates.generation", -0.1, 0.1);
  
  std::vector<SimpleCut> lambdacuts
  {
    EqualsCut("Candidates.pid",              3122               ),
//     RangeCut ("Candidates.chi2_prim_first",  18.4207, HugeValue ),
//     RangeCut ("Candidates.chi2_prim_second", 18.4207, HugeValue ),
//     RangeCut ("Candidates.distance",         0.,      1.        ),
//     RangeCut ("Candidates.chi2_geo",         0.,      3.        ),
//     RangeCut ("Candidates.l_over_dl",        0.,      5.        ),
//     RangeCut ("Candidates.cosine_first",     0.,      HugeValue ),
//     RangeCut ("Candidates.cosine_second",    0.,      HugeValue ),
//     RangeCut ("Candidates.chi2_topo",        0.,      15.       ),
//     RangeCut ("Candidates.chi2prim_mother",  0.,      15.       ),
  };
  
  std::vector<SimpleCut> xicuts
  {
    EqualsCut("Candidates.pid",              3312               ),
//     RangeCut ("Candidates.chi2_prim_first",  18.4207, HugeValue ),
//     RangeCut ("Candidates.chi2_prim_second", 18.4207, HugeValue ),
//     RangeCut ("Candidates.distance",         0.,      1.        ),
//     RangeCut ("Candidates.chi2_geo",         0.,      3.        ),
//     RangeCut ("Candidates.l_over_dl",        0.,      5.        ),
//     RangeCut ("Candidates.cosine_first",     0.,      HugeValue ),
//     RangeCut ("Candidates.cosine_second",    0.,      HugeValue ),
//     RangeCut ("Candidates.chi2_topo",        0.,      15.       ),
//     RangeCut ("Candidates.chi2prim_mother",  0.,      15.       ),
  };
  
  struct QaSetup
  {
    std::string name_;
    SimpleCut cut_;
  };
  
  std::vector<QaSetup> qasetups
  {
    {"All", allcut},
    {"Sgnl", sgnlcut},
    {"Bckgr", bckgrcut}
  };
  
  struct QaPid
  {
    std::string name_;
    std::vector<SimpleCut> cuts_;
    std::string axisname_;
    std::string daughter_name_1_;
    std::string daughter_name_2_;
  };
  
  std::vector<QaPid> qapids
  {
//     {"Lambda", lambdacuts, "#Lambda", "#pi^{-}", "p"},
    {"Ksi",    xicuts,     "#Xi^{-}", "#pi^{-}", "#Lambda"}
  };
  
  for(auto& qp : qapids)
    for(auto& qs : qasetups)
    {
      auto current_cuts = qp.cuts_;
      current_cuts.push_back(qs.cut_);
      
      Cuts* selection_cuts = new Cuts((qs.name_ + "_" + qp.name_).c_str(), current_cuts);
      float mass_axis_up;
      if(qp.name_ == "Lambda")
        mass_axis_up = 2.;
      else
        mass_axis_up = 3.;
      
      //**************************************** 1D histograms of Candidates ******************************************************************************************************************************************************
      task->AddH1({("m_{" + qp.axisname_ + "}, GeV/c^{2}").c_str(),                                Variable::FromString("Candidates.mass"),             {nbins, 1, mass_axis_up}}, selection_cuts);
      task->AddH1({("p_{" + qp.axisname_ + "}, GeV/c").c_str(),                                    Variable::FromString("Candidates.p"),                {100, 0,   20}},           selection_cuts);
      task->AddH1({("p_{X " + qp.axisname_ + "}, GeV/c").c_str(),                                  Variable::FromString("Candidates.px"),               {100, -5,   5}},           selection_cuts);
      task->AddH1({("p_{Y " + qp.axisname_ + "}, GeV/c").c_str(),                                  Variable::FromString("Candidates.py"),               {100, -5,   5}},           selection_cuts);
      task->AddH1({("p_{Z " + qp.axisname_ + "}, GeV/c").c_str(),                                  Variable::FromString("Candidates.pz"),               {100, 0,   20}},           selection_cuts);
      task->AddH1({("p_{T " + qp.axisname_ + "}, GeV/c").c_str(),                                  Variable::FromString("Candidates.pT"),               {100, 0,   5 }},           selection_cuts);
      task->AddH1({("y_{LAB " + qp.axisname_ + "}").c_str(),                                       Variable::FromString("Candidates.rapidity"),         {40,  0,   4 }},           selection_cuts);
      task->AddH1({("#varphi_{" + qp.axisname_ + "}, rad").c_str(),                                Variable::FromString("Candidates.phi"),              {100, -PI, PI}},           selection_cuts);
      task->AddH1({("X_{" + qp.axisname_ + "}, cm").c_str(),                                       Variable::FromString("Candidates.x"),                {200, -50, 50}},           selection_cuts);
      task->AddH1({("Y_{" + qp.axisname_ + "}, cm").c_str(),                                       Variable::FromString("Candidates.y"),                {200, -50, 50}},           selection_cuts);
      task->AddH1({("Z_{" + qp.axisname_ + "}, cm").c_str(),                                       Variable::FromString("Candidates.z"),                {360, -10, 80}},           selection_cuts);
      task->AddH1({"generation",                                                                   Variable::FromString("Candidates.generation"),       {5, -1.5, 3.5}},           selection_cuts);
      task->AddH1({("#chi^{2}_{prim, " + qp.daughter_name_1_ + "}").c_str(),                       Variable::FromString("Candidates.chi2_prim_first"),  {nbins, 0, 400}},          selection_cuts);
      task->AddH1({("#chi^{2}_{prim, " + qp.daughter_name_2_ + "}").c_str(),                       Variable::FromString("Candidates.chi2_prim_second"), {nbins, 0, 400}},          selection_cuts);
      task->AddH1({"DCA, cm",                                                                      Variable::FromString("Candidates.distance"),         {nbins, 0, 20}},           selection_cuts);
      task->AddH1({("cos(#alpha_{" + qp.axisname_ + ", " + qp.daughter_name_1_ + "})").c_str(),    Variable::FromString("Candidates.cosine_first"),     {nbins, 0.5, 1}},          selection_cuts);
      task->AddH1({("cos(#alpha_{" + qp.axisname_ + ", " + qp.daughter_name_2_ + "})").c_str(),    Variable::FromString("Candidates.cosine_second"),    {nbins, 0.5, 1}},          selection_cuts);
      task->AddH1({"#chi^{2}_{geo}",                                                               Variable::FromString("Candidates.chi2_geo"),         {nbins, 0, 100}},          selection_cuts);
      task->AddH1({"L/#Delta L",                                                                   Variable::FromString("Candidates.l_over_dl"),        {nbins, 0, 100}},          selection_cuts);
      task->AddH1({"#chi^{2}_{topo}",                                                              Variable::FromString("Candidates.chi2_topo"),        {nbins, 0, 100}},          selection_cuts);
      task->AddH1({("#chi^{2}_{prim, " + qp.axisname_ + "}").c_str(),                              Variable::FromString("Candidates.chi2_prim_mother"), {nbins, 0, 100}},          selection_cuts);
      task->AddH1({("#Delta m_{" + qp.axisname_ + "} / #sigma_{m, " + qp.axisname_ + "}").c_str(), Variable::FromString("Candidates.invmass_discr"),    {nbins, 0, 100}},          selection_cuts);
      task->AddH2({("#chi^{2}_{prim, " + qp.axisname_ + "}").c_str(), Variable::FromString("Candidates.chi2_prim_mother"), {nbins, 0, 100}}, {"#chi^{2}_{topo}", Variable::FromString("Candidates.chi2_topo"), {nbins, 0, 100}}, selection_cuts);
      //***************************************************************************************************************************************************************************************************************************
      
      if(qs.name_ != "Sgnl") continue;
      
      //***************************************** 2D histograms Candidates-Simulated **********************************************************************************************************************************************
      task->AddH2({("p^{sim}_{" + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Simulated.p"), {100, 0, 20}}, {("p^{reco}_{" + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Candidates.p"), {100, 0, 20}}, selection_cuts);
      task->AddH2({("p^{sim}_{T " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Simulated.pT"), {100, 0, 5}}, {("p^{reco}_{T " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Candidates.pT"), {100, 0, 5}}, selection_cuts);
      task->AddH2({("p^{sim}_{X " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Simulated.px"), {100, -5, 5}}, {("p^{reco}_{X " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Candidates.px"), {100, -5, 5}}, selection_cuts);
      task->AddH2({("p^{sim}_{Y " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Simulated.py"), {100, -5, 5}}, {("p^{reco}_{Y " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Candidates.py"), {100, -5, 5}}, selection_cuts);
      task->AddH2({("p^{sim}_{Z " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Simulated.pz"), {100, 0, 20}}, {("p^{reco}_{Z " + qp.axisname_ + "}, GeV/c").c_str(), Variable::FromString("Candidates.pz"), {100, 0, 20}}, selection_cuts);
      task->AddH2({("y^{sim}_{LAB " + qp.axisname_ + "}").c_str(), Variable::FromString("Simulated.rapidity"), {40, 0, 4}}, {("y^{reco}_{LAB " + qp.axisname_ + "}").c_str(), Variable::FromString("Candidates.rapidity"), {40, 0, 4}}, selection_cuts);
      task->AddH2({("#varphi^{sim}_{" + qp.axisname_ + "}").c_str(), Variable::FromString("Simulated.phi"), {100, -PI, PI}}, {("#varphi^{reco}_{" + qp.axisname_ + "}").c_str(), Variable::FromString("Candidates.phi"), {100, -PI, PI}}, selection_cuts);
      task->AddH2({("X^{sim}_{" + qp.axisname_ + "}, cm").c_str(), Variable::FromString("Simulated.x"), {200, -50, 50}}, {("X^{reco}_{" + qp.axisname_ + "}, cm").c_str(), Variable::FromString("Candidates.x"), {200, -50, 50}}, selection_cuts);
      task->AddH2({("Y^{sim}_{" + qp.axisname_ + "}, cm").c_str(), Variable::FromString("Simulated.y"), {200, -50, 50}}, {("Y^{reco}_{" + qp.axisname_ + "}, cm").c_str(), Variable::FromString("Candidates.y"), {200, -50, 50}}, selection_cuts);
      task->AddH2({("Z^{sim}_{" + qp.axisname_ + "}, cm").c_str(), Variable::FromString("Simulated.z"), {360, -10, 80}}, {("Z^{reco}_{" + qp.axisname_ + "}, cm").c_str(), Variable::FromString("Candidates.z"), {360, -10, 80}}, selection_cuts);    
      //***************************************************************************************************************************************************************************************************************************
      
      //***************************************** 1D histograms difference Candidates-Simulated ***********************************************************************************************************************************
      Variable diff_p("diff_p", {{"Candidates", "p"}, {"Simulated", "p"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("p^{reco}_{" + qp.axisname_ + "} - p^{sim}_{" + qp.axisname_ + "}, GeV/c").c_str(), diff_p, {100, -1, 1}}, selection_cuts);
      Variable diff_pX("diff_px", {{"Candidates", "px"}, {"Simulated", "px"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("p^{reco}_{X " + qp.axisname_ + "} - p^{sim}_{X " + qp.axisname_ + "}, GeV/c").c_str(), diff_pX, {100, -0.1, 0.1}}, selection_cuts);
      Variable diff_pY("diff_py", {{"Candidates", "py"}, {"Simulated", "py"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("p^{reco}_{Y " + qp.axisname_ + "} - p^{sim}_{Y " + qp.axisname_ + "}, GeV/c").c_str(), diff_pY, {100, -0.1, 0.1}}, selection_cuts);
      Variable diff_pZ("diff_pz", {{"Candidates", "pz"}, {"Simulated", "pz"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("p^{reco}_{Z " + qp.axisname_ + "} - p^{sim}_{Z " + qp.axisname_ + "}, GeV/c").c_str(), diff_pZ, {100, -0.5, 0.5}}, selection_cuts);
      Variable diff_pT("diff_pT", {{"Candidates", "pT"}, {"Simulated", "pT"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("p^{reco}_{T " + qp.axisname_ + "} - p^{sim}_{T " + qp.axisname_ + "}, GeV/c").c_str(), diff_pT, {100, -0.1, 0.1}}, selection_cuts);
      Variable diff_rapidity("diff_rapidity", {{"Candidates", "rapidity"}, {"Simulated", "rapidity"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("y^{reco}_{LAB " + qp.axisname_ + "} - y^{sim}_{LAB " + qp.axisname_ + "}").c_str(), diff_rapidity, {100, -0.1, 0.1}}, selection_cuts);
      Variable diff_phi("diff_phi", {{"Candidates", "phi"}, {"Simulated", "phi"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("#varphi^{reco}_{" + qp.axisname_ + "} - #varphi^{sim}_{" + qp.axisname_ + "}, rad").c_str(), diff_phi, {100, -0.1, 0.1}}, selection_cuts);
      Variable diff_x("diff_x", {{"Candidates", "x"}, {"Simulated", "x"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("X^{reco}_{" + qp.axisname_ + "} - X^{sim}_{" + qp.axisname_ + "}, cm").c_str(), diff_x, {200, -20, 20}}, selection_cuts);
      Variable diff_y("diff_y", {{"Candidates", "y"}, {"Simulated", "y"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("Y^{reco}_{" + qp.axisname_ + "} - Y^{sim}_{" + qp.axisname_ + "}, cm").c_str(), diff_y, {200, -20, 20}}, selection_cuts);
      Variable diff_z("diff_z", {{"Candidates", "z"}, {"Simulated", "z"}}, []( std::vector<double>& var ) { return var.at(0)-var.at(1); });
      task->AddH1({("Z^{reco}_{" + qp.axisname_ + "} - Z^{sim}_{" + qp.axisname_ + "}, cm").c_str(), diff_z, {360, -10, 80}}, selection_cuts);
      //***************************************************************************************************************************************************************************************************************************
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