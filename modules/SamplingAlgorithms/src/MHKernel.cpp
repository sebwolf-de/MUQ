#include "MUQ/SamplingAlgorithms/MHKernel.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MHKernel)

MHKernel::MHKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : MCMCKernel(pt, problem) {}

MHKernel::~MHKernel() {}

void MHKernel::EvaluateImpl(ref_vector<boost::any> const& inputs) {

  boost::any temp = inputs.at(0);

  // current state
  std::shared_ptr<SamplingState> current = boost::any_cast<std::shared_ptr<SamplingState> >(inputs.at(0).get());
  assert(current);

  // propose a new point
  std::shared_ptr<SamplingState> prop = boost::any_cast<std::shared_ptr<SamplingState> >(proposal->Sample(current));

  // compute acceptance probability
  if( !current->HasMeta("LogTarget") ){
     current->meta["LogTarget"] = problem->EvaluateLogTarget(ref_vector<boost::any>(1, std::cref(current->state.at(0))));
  }
  if( !prop->HasMeta("LogTarget") ){
     prop->meta["LogTarget"] = problem->EvaluateLogTarget(ref_vector<boost::any>(1, std::cref(prop->state.at(0))));
  }

  double propTarget = AnyCast(prop->meta["LogTarget"]);
  double currentTarget = AnyCast(current->meta["LogTarget"]);

  const double alpha = std::exp(propTarget + proposal->LogDensity(current, prop)- currentTarget -proposal->LogDensity(prop, current));

  // accept/reject
  outputs.resize(1);
  if( RandomGenerator::GetUniform()<alpha ) {
    outputs[0] = prop;
  } else {
    current->weight += 1.0;
    outputs[0] = boost::none;
  }
}
