#include "CutFactory.hpp"

#include <stdexcept>

namespace AnalysisTree {
namespace QA {

CutFactory& CutFactory::Instance() {
  static CutFactory instance;
  return instance;
}

void CutFactory::SetSharedCutsNode(const YAML::Node& shared_cuts_node) {
  // Assigning an undefined node into an already-constructed Node member via operator= throws yaml-cpp's InvalidNode (unlike fresh-initializing a local one): only assign when there's actually something defined; the
  // member's own default constructor already leaves it equivalently "undefined" otherwise.
  if (shared_cuts_node.IsDefined()) {
    shared_cuts_node_ = shared_cuts_node;
  }
}

void CutFactory::RegisterCutType(const std::string& type_name, CutCreatorFunction creator) {
  const auto result = cut_creators_.emplace(type_name, std::move(creator));
  if (!result.second) {
    throw std::runtime_error("CutFactory::RegisterCutType(): custom cut type '" + type_name + "' is already registered");
  }
}

AnalysisTree::SimpleCut CutFactory::CreateCustomCut(const std::string& type_name, const YAML::Node& params) const {
  const auto it = cut_creators_.find(type_name);
  if (it == cut_creators_.end()) {
    throw std::runtime_error("CutFactory::CreateCustomCut(): unknown custom cut type '" + type_name
                             + "' (is its .cpp file linked into this binary? see services/README.md)");
  }
  return it->second(params);
}

AnalysisTree::SimpleCut CutFactory::BuildSimpleCut(const YAML::Node& simple_cut_node) const {
  const auto type = simple_cut_node["type"].as<std::string>();
  const auto title = GetOrDefault<std::string>(simple_cut_node["title"], "");

  if (type == "range") {
    return RangeCut(simple_cut_node["variable"].as<std::string>(),
                    simple_cut_node["min"].as<double>(),
                    simple_cut_node["max"].as<double>(),
                    title);
  }
  if (type == "equals") {
    return EqualsCut(simple_cut_node["variable"].as<std::string>(), simple_cut_node["value"].as<int>(), title);
  }
  if (type == "not_equals") {
    const auto variable = simple_cut_node["variable"].as<std::string>();
    const int value = simple_cut_node["value"].as<int>();
    return AnalysisTree::SimpleCut(
        {variable}, [value](std::vector<double> v) { return static_cast<int>(v[0]) != value; }, title);
  }
  if (type == "custom") {
    return CreateCustomCut(simple_cut_node["name"].as<std::string>(), simple_cut_node);
  }
  if (type == "or") {
    return BuildOrCut(simple_cut_node);
  }
  throw std::runtime_error("CutFactory::BuildSimpleCut(): unknown simple cut type '" + type
                           + "' (expected 'range', 'equals', 'not_equals', 'custom' or 'or')");
}

AnalysisTree::SimpleCut CutFactory::BuildOrCut(const YAML::Node& or_node) const {
  const auto title = GetOrDefault<std::string>(or_node["title"], "");
  const YAML::Node clauses = or_node["clauses"];
  if (!clauses.IsDefined() || !clauses.IsSequence() || clauses.size() == 0) {
    throw std::runtime_error("CutFactory::BuildOrCut(): 'or' requires a non-empty 'clauses' list");
  }

  // Each clause references exactly one variable, so clause predicate i is evaluated against v[i] below: order matters, hence keeping "variables" and "predicates" as two parallel vectors instead of e.g. a map.
  std::vector<std::string> variables;
  std::vector<std::function<bool(double)>> predicates;
  variables.reserve(clauses.size());
  predicates.reserve(clauses.size());

  for (const auto& clause : clauses) {
    const auto clause_type = clause["type"].as<std::string>();
    variables.push_back(clause["variable"].as<std::string>());

    if (clause_type == "range") {
      const auto min = clause["min"].as<double>();
      const auto max = clause["max"].as<double>();
      predicates.push_back([min, max](double x) { return x >= min && x <= max; });
    } else if (clause_type == "equals") {
      const int value = clause["value"].as<int>();
      predicates.push_back([value](double x) { return static_cast<int>(x) == value; });
    } else if (clause_type == "not_equals") {
      const int value = clause["value"].as<int>();
      predicates.push_back([value](double x) { return static_cast<int>(x) != value; });
    } else {
      throw std::runtime_error("CutFactory::BuildOrCut(): unsupported clause type '" + clause_type
                               + "' inside 'or' (expected 'range', 'equals' or 'not_equals')");
    }
  }

  return AnalysisTree::SimpleCut(
      variables,
      [predicates](std::vector<double> v) {
        for (std::size_t i = 0; i < predicates.size(); ++i) {
          if (predicates[i](v[i])) return true;
        }
        return false;
      },
      title);
}

AnalysisTree::Cuts* CutFactory::BuildCuts(const YAML::Node& cuts_field, const std::string& default_name) const {
  if (!cuts_field.IsDefined() || cuts_field.IsNull()) {
    return nullptr;
  }

  if (cuts_field.IsScalar()) {
    const auto name = cuts_field.as<std::string>();
    const YAML::Node shared_entry = shared_cuts_node_[name];
    if (!shared_entry.IsDefined()) {
      throw std::runtime_error("CutFactory::BuildCuts(): no shared cut named '" + name
                               + "' found under the top-level 'cuts:' map");
    }
    return BuildCuts(shared_entry, name);
  }

  const YAML::Node simple_cuts = cuts_field["simple_cuts"];
  if (!simple_cuts.IsDefined() || !simple_cuts.IsSequence()) {
    throw std::runtime_error("CutFactory::BuildCuts(): expected a 'simple_cuts' list while building cut '"
                             + default_name + "'");
  }

  std::vector<AnalysisTree::SimpleCut> cuts;
  cuts.reserve(simple_cuts.size());
  for (const auto& simple_cut_node : simple_cuts) {
    cuts.emplace_back(BuildSimpleCut(simple_cut_node));
  }
  return new AnalysisTree::Cuts(default_name, cuts);
}

}// namespace QA
}// namespace AnalysisTree
