#include "MUQ/SamplingAlgorithms/MCMCKernel.h"

using namespace muq::SamplingAlgorithms;

MCMCKernel::MCMCKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem), proposal(MCMCProposal::Construct(pt)) {}

MCMCKernel::~MCMCKernel() {}

void MCMCKernel::PostStep(unsigned int const t, std::shared_ptr<SamplingState> state) {
  proposal->Adapt(t, state);
}
