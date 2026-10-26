#ifndef ANALYSISTREEQA_SERVICES_STEER_CUTFACTORY_HPP_
#define ANALYSISTREEQA_SERVICES_STEER_CUTFACTORY_HPP_

#include <functional>
#include <map>
#include <string>

#include <yaml-cpp/yaml.h>

#include "AnalysisTree/Cuts.hpp"
#include "AnalysisTree/SimpleCut.hpp"

namespace AnalysisTree {
namespace QA {

/// Like YAML::Node::as<T>(default_value), but never touches yaml-cpp's own conversion machinery for a node that is entirely absent (as opposed to present with an unexpected type,
/// which should still surface as an error). Every task/cut YAML field that's documented as "optional" should go through this rather than the raw .as<T>(default) overload.
template<typename T>
T GetOrDefault(const YAML::Node& node, const T& default_value) {
  return node.IsDefined() ? node.as<T>() : default_value;
}

/// Builds AnalysisTree::Cuts from YAML, either declaratively (range/equals/not_equals, optionally combined with logical OR via `or: {clauses: [...]}`) or via a named,
/// C++-registered custom predicate (for logic that cannot be expressed declaratively at all, e.g. a ratio of two fields), see services/README.md for the full schema and examples.
/// Every entry in one `simple_cuts:` list (including an `or` entry) is always combined with the others via logical AND, that's the one combination AnalysisTree::Cuts itself provides.
class CutFactory {
 public:
  /// Builds one AnalysisTree::SimpleCut from a "custom" cut's YAML params.
  using CutCreatorFunction = std::function<AnalysisTree::SimpleCut(const YAML::Node& params)>;

  static CutFactory& Instance();

  /// Called once by Run at startup with the top-level "cuts:" map, so that a scalar-string cut reference can be resolved wherever it's encountered.
  void SetSharedCutsNode(const YAML::Node& shared_cuts_node);

  /// cuts_field is either: a scalar string (a name looked up in the shared cuts map set via SetSharedCutsNode), an inline mapping ({simple_cuts: [...]}), or a missing/null node (-> nullptr).
  /// This is what every task should call for its own "cuts:" YAML field. The returned Cuts* is heap-allocated and owned by the caller (matching how AnalysisTree::Cuts* is used throughout the rest of AnalysisTreeQA).
  AnalysisTree::Cuts* BuildCuts(const YAML::Node& cuts_field, const std::string& default_name) const;

  void RegisterCutType(const std::string& type_name, CutCreatorFunction creator);

 private:
  CutFactory() = default;

  AnalysisTree::SimpleCut BuildSimpleCut(const YAML::Node& simple_cut_node) const;
  AnalysisTree::SimpleCut CreateCustomCut(const std::string& type_name, const YAML::Node& params) const;
  /// {type: or, clauses: [...]}: combines several range/equals/not_equals clauses (each on its own variable) with logical OR into one SimpleCut. AnalysisTree::Cuts itself only ever ANDs its SimpleCuts, so this is the
  /// only way to express OR: as a single, self-contained SimpleCut.
  AnalysisTree::SimpleCut BuildOrCut(const YAML::Node& or_node) const;

  YAML::Node shared_cuts_node_;
  std::map<std::string, CutCreatorFunction> cut_creators_;
};

}// namespace QA
}// namespace AnalysisTree

/// Registers a named, C++-defined cut predicate once, at static-init time. Lambda must have signature: AnalysisTree::SimpleCut(const YAML::Node&). Referenced from YAML as: {type: custom, name: TypeName, ...params}.
#define ATQA_REGISTER_CUT(TypeName, Lambda)                                                                      \
  namespace {                                                                                                    \
  struct TypeName##_CutRegistrar {                                                                               \
    TypeName##_CutRegistrar() { ::AnalysisTree::QA::CutFactory::Instance().RegisterCutType(#TypeName, Lambda); } \
  } TypeName##_cut_registrar_instance;                                                                           \
  }

#endif//ANALYSISTREEQA_SERVICES_STEER_CUTFACTORY_HPP_
