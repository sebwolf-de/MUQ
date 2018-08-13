#include "MUQ/SamplingAlgorithms/GMHKernel.h"

#include <parcer/Queue.h>

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(GMHKernel)

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : MHKernel(pt, problem), N(pt.get<unsigned int>("NumProposals")) {}

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) : MHKernel(pt, problem, proposalIn), N(pt.get<unsigned int>("NumProposals")) {}

GMHKernel::~GMHKernel() {}

void GMHKernel::PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) {
  // propose N steps
  std::cout << "propose N steps" << std::endl;
  std::vector<std::shared_ptr<SamplingState> > proposals(N+1);
  proposals[0] = state;
  for( auto it = proposals.begin()+1; it!=proposals.end(); ++it ) {
    *it = proposal->Sample(state); 
  }

  // compute log-target
  std::cout << "compute log-target" << std::endl;
  for( auto it : proposals ) { std::cout << problem->LogDensity(it) << std::endl; }

  // compute stationary acceptance transition probability
  std::cout << "compute acceptance ratio" << std::endl;

}

/*std::vector<std::shared_ptr<SamplingState> > GMHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> state) {
  return std::vector<std::shared_ptr<SamplingState> >(1, nullptr);
  }*/

