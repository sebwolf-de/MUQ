#include "MUQ/SamplingAlgorithms/GMHKernel.h"

#include <parcer/Queue.h>

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(GMHKernel)

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : MHKernel(pt, problem) {}

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) : MHKernel(pt, problem, proposalIn) {}

GMHKernel::~GMHKernel() {}

void GMHKernel::PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) {}

std::vector<std::shared_ptr<SamplingState> > GMHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> state) {
  return std::vector<std::shared_ptr<SamplingState> >(1, nullptr);
}

