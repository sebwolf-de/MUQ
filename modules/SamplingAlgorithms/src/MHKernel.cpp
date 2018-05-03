#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MHKernel)

MHKernel::MHKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem),
                                                                                    blockInd(pt.get("BlockIndex",0)),
                                                                                    proposal(MCMCProposal::Construct(pt, problem)) {}


MHKernel::~MHKernel() {}

std::shared_ptr<SamplingState> MHKernel::Step(std::shared_ptr<SamplingState> prevState){

  assert(proposal);

  // propose a new point
  boost::any propAny = proposal->Sample(prevState);
  std::shared_ptr<SamplingState> prop = AnyCast(propAny);

  // compute acceptance probability
  double propTarget;
  double currentTarget;

  if( prevState->HasMeta("LogTarget") ){
    currentTarget = AnyCast( prevState->meta["LogTarget"]);
  }else{
    currentTarget = problem->LogDensity(prevState);
    prevState->meta["LogTarget"] = currentTarget;
  }

  propTarget = problem->LogDensity(prop);
  prop->meta["LogTarget"] = propTarget;

  // Aceptance probability
  const double forwardPropDens = proposal->LogDensity(prevState, prop);
  const double backPropDens = proposal->LogDensity(prop, prevState);
  const double alpha = std::exp(propTarget + forwardPropDens - currentTarget - backPropDens);

  // accept/reject
  if( RandomGenerator::GetUniform()<alpha ) {
    return prop;
  } else {
    prevState->weight += 1.0;
    return prevState;
  }
}
