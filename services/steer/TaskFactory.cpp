#include "TaskFactory.hpp"

#include <stdexcept>

namespace AnalysisTree {
namespace QA {

TaskFactory& TaskFactory::Instance() {
  static TaskFactory instance;
  return instance;
}

void TaskFactory::RegisterTaskType(const std::string& type_name, TaskCreatorFunction creator) {
  const auto result = task_creators_.emplace(type_name, std::move(creator));
  if (!result.second) {
    throw std::runtime_error("TaskFactory::RegisterTaskType(): task type '" + type_name + "' is already registered");
  }
}

AnalysisTree::Task* TaskFactory::CreateTask(const std::string& type_name, const YAML::Node& node) const {
  const auto it = task_creators_.find(type_name);
  if (it == task_creators_.end()) {
    throw std::runtime_error("TaskFactory::CreateTask(): unknown task type '" + type_name
                             + "' (is its .cpp file linked into this binary? see services/README.md)");
  }
  return it->second(node);
}

std::vector<std::string> TaskFactory::RegisteredTypeNames() const {
  // task_creators_ is a std::map<std::string, ...>, which already iterates
  // in sorted key order - no extra sort needed.
  std::vector<std::string> names;
  names.reserve(task_creators_.size());
  for (const auto& entry : task_creators_) {
    names.push_back(entry.first);
  }
  return names;
}

}// namespace QA
}// namespace AnalysisTree
