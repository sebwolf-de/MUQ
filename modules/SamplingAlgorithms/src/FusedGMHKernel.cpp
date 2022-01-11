#include "MUQ/SamplingAlgorithms/FusedGMHKernel.h"

#include <Eigen/Eigenvalues>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/AnyHelpers.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(FusedGMHKernel)

FusedGMHKernel::FusedGMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : GMHKernel(pt, problem) {}

FusedGMHKernel::FusedGMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) :
  GMHKernel(pt, problem, proposalIn) {}

FusedGMHKernel::~FusedGMHKernel() {}

void GMHKernel::SerialProposal(unsigned int const t, std::shared_ptr<SamplingState> state) {

  std::shared_ptr<SamplingState> helpState;
  helpState->state.resize(N, nullptr); // TODO: check for correct pointer syntax

  // If the current state does not have LogTarget information, add it
  if(! state->HasMeta("LogTarget"))
    state->meta["LogTarget"] = 0.0; // dummy value to avoid unfuded sim ... problem->LogDensity(state);


  // propose the points
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;

  for(unsigned int j = 0; j<N; j++) {
    helpState->state.at(j) = proposal->Sample(state); // TODO: check for correct pointer syntax
  }

  // Run fused simulation
  problem->LogDensity(helpState);

  // Transfer LogDensity data to proposedStates
  unsigned int j = 0;
  for(auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) {
    *it = helpState->state.at(j++);
    (*it)->meta["LogTarget"] = helpState->meta["LogTarget"][j++];
  }

  // evaluate the target density
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  for( unsigned int i=0; i<Np1; ++i )
    R(i) = AnyCast(proposedStates[i]->meta["LogTarget"]);

  // compute stationary transition probability
  AcceptanceDensity(R);
}