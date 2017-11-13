#include "MUQ/SamplingAlgorithms/MCMC.h"

using namespace muq::SamplingAlgorithms;

MCMC::MCMC() : SamplingAlgorithm() {}

MCMC::~MCMC() {}

std::shared_ptr<TransitionKernel> MCMC::Kernel(boost::property_tree::ptree& pt, std::shared_ptr<SamplingProblem> problem) const {
  return TransitionKernel::Construct(pt, problem);
}
