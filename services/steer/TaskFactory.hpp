#ifndef ANALYSISTREEQA_SERVICES_STEER_TASKFACTORY_HPP_
#define ANALYSISTREEQA_SERVICES_STEER_TASKFACTORY_HPP_

#include <functional>
#include <map>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "AnalysisTree/Task.hpp"

namespace AnalysisTree {
namespace QA {

/// Creates concrete AnalysisTree::Task instances by a YAML-supplied "task" name. Task types (built-in ones under services/steer/tasks/, or private ones a user drops into the top-level tasks/ directory) register
/// themselves once via ATQA_REGISTER_TASK - no central list needs editing to add a new task type.
class TaskFactory {
 public:
  /// Builds one AnalysisTree::Task from its own YAML node (the full entry under the top-level "tasks:" list, including "task"/"name"/"cuts"/etc.).
  using TaskCreatorFunction = std::function<AnalysisTree::Task*(const YAML::Node&)>;

  static TaskFactory& Instance();

  void RegisterTaskType(const std::string& type_name, TaskCreatorFunction creator);
  AnalysisTree::Task* CreateTask(const std::string& type_name, const YAML::Node& node) const;

  /// Every task type name currently registered (i.e. whose .cpp file got linked into this binary), in sorted order for `--print-registered-tasks`.
  std::vector<std::string> RegisteredTypeNames() const;

 private:
  TaskFactory() = default;
  std::map<std::string, TaskCreatorFunction> task_creators_;
};

/// Registers ConcreteTask under type_name at static-initialization time. Instantiated once per task type via ATQA_REGISTER_TASK below.
template<typename ConcreteTask>
struct TaskTypeRegistrar {
  explicit TaskTypeRegistrar(const std::string& type_name) {
    TaskFactory::Instance().RegisterTaskType(
        type_name, [](const YAML::Node& node) -> AnalysisTree::Task* { return new ConcreteTask(node); });
  }
};

}// namespace QA
}// namespace AnalysisTree

/// Place once in a concrete task's .cpp file, INSIDE the same namespace the class is declared in, using its unqualified name (ClassName must have a `explicit ClassName(const YAML::Node&)` constructor), e.g. for a class
/// AnalysisTree::QA::TrackQaTask, call ATQA_REGISTER_TASK(TrackQaTask) from within `namespace AnalysisTree { namespace QA { ... } }`, not after it with the qualified name (breaks the token-pasted registrar variable name
/// and would register the wrong, fully-qualified type string). Linking that .cpp's object file into the final binary is sufficient for ClassName to become selectable from YAML via `task: ClassName`, see services/README.md.
#define ATQA_REGISTER_TASK(ClassName)                                                             \
  namespace {                                                                                     \
  static ::AnalysisTree::QA::TaskTypeRegistrar<ClassName> ClassName##_task_registrar(#ClassName); \
  }

#endif//ANALYSISTREEQA_SERVICES_STEER_TASKFACTORY_HPP_
