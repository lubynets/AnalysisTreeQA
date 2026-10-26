#include "Run.hpp"

#include <filesystem>
#include <fstream>
#include <unistd.h>
#include <utility>

#include <TFile.h>

#include "AnalysisTree/TaskManager.hpp"

#include "Task.hpp"

#include "CutFactory.hpp"
#include "TaskFactory.hpp"

namespace AnalysisTree {
namespace QA {

namespace {

// AnalysisTree::TaskManager::Init()/AnalysisTree::Chain::InitChain() treats every input path as a "filelist" text file (one ROOT file path per line): passing a ROOT file directly there makes it try to parse binary content as
// text, which fails badly (see services/README.md). Detect that case by the ROOT file magic bytes ("root", at the very start of the file) rather than the extension, so both forms work interchangeably via -i/input.files.
bool LooksLikeRootFile(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return false;
  }
  char magic[4] = {};
  in.read(magic, sizeof(magic));
  return in.gcount() == static_cast<std::streamsize>(sizeof(magic)) && magic[0] == 'r' && magic[1] == 'o'
      && magic[2] == 'o' && magic[3] == 't';
}

// Transparently wraps a direct ROOT file path in a one-line temporary filelist, so callers can pass either a ROOT file or an actual filelist via -i/input.files. Returns `path` unchanged if it isn't a ROOT file (already
// assumed to be a filelist). Appends any created temp file to `temp_files_to_clean_up` so the caller can remove it once the run is done.
std::string ResolveAsFilelist(const std::string& path, std::vector<std::string>& temp_files_to_clean_up) {
  if (!LooksLikeRootFile(path)) {
    return path;
  }

  const auto temp_path = (std::filesystem::temp_directory_path()
                          / ("analysistreeqa_filelist_" + std::to_string(getpid()) + "_"
                             + std::to_string(temp_files_to_clean_up.size()) + ".txt"))
                             .string();
  std::ofstream out(temp_path);
  out << std::filesystem::absolute(path).string() << "\n";
  temp_files_to_clean_up.push_back(temp_path);
  return temp_path;
}

}// namespace

Run::Run(Config config) : config_(std::move(config)) {}

void Run::Exec() {
  auto* out_file = new TFile(config_.OutputFile().c_str(), config_.Overwrite() ? "recreate" : "create");

  CutFactory::Instance().SetSharedCutsNode(config_.SharedCutsNode());

  auto* man = AnalysisTree::TaskManager::GetInstance();

  for (const auto& entry : config_.TaskNodes()) {
    const auto type_name = entry["task"].as<std::string>();
    const YAML::Node name_node = entry["name"];
    const auto task_name = name_node.IsDefined() ? name_node.as<std::string>() : type_name;

    auto* task = TaskFactory::Instance().CreateTask(type_name, entry);

    if (auto* qa_task = dynamic_cast<AnalysisTree::QA::Task*>(task)) {
      qa_task->AttachOutputFile(out_file);
      qa_task->SetTopLevelDirName(task_name);
    }

    man->AddTask(task);
  }

  std::vector<std::string> temp_filelists_to_clean_up;
  std::vector<std::string> effective_input_files;
  effective_input_files.reserve(config_.InputFiles().size());
  for (const auto& input_file : config_.InputFiles()) {
    effective_input_files.push_back(ResolveAsFilelist(input_file, temp_filelists_to_clean_up));
  }

  man->Init(effective_input_files, {config_.TreeName()});
  man->Run(config_.NEvents());
  man->Finish();

  out_file->Write();
  out_file->Close();

  for (const auto& temp_file : temp_filelists_to_clean_up) {
    std::filesystem::remove(temp_file);
  }
}

}// namespace QA
}// namespace AnalysisTree
