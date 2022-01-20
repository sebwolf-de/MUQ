#include "MUQ/SamplingAlgorithms/FusedGMHKernel.h"

#include <Eigen/Eigenvalues>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/AnyHelpers.h"

#include <iostream>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(FusedGMHKernel)

FusedGMHKernel::FusedGMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : 
  GMHKernel(pt, problem) {}

FusedGMHKernel::FusedGMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) :
  GMHKernel(pt, problem, proposalIn) {}

FusedGMHKernel::~FusedGMHKernel() {}

void FusedGMHKernel::PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) {
  // propose N steps
  FusedProposal(t, state);
}

void FusedGMHKernel::FusedProposal(unsigned int const t, std::shared_ptr<SamplingState> state) {

  std::shared_ptr<SamplingState> helpState;
  helpState->state.resize(0); // TODO: check for correct pointer syntax

  // If the current state does not have LogTarget information, add it
  if(! state->HasMeta("LogTarget"))
    state->meta["LogTarget"] = 0.0; // dummy value to avoid unfuded sim ... problem->LogDensity(state);


  // propose the points
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;

  for(unsigned int j = 0; j<N; j++) {
    helpState->state.push_back(proposal->Sample(state)->state[0]);
  }
  std::cout << "Size of N: " << N << std::endl;
  std::cout << "Size of state vector: " << helpState->state.size() << std::endl;
  // for(unsigned int j = 0; j<N; j++) {
  //   helpState->state.at(j) = proposal->Sample(state)->state[0];
  // }
  // for(auto it = helpState->state.begin()+1; it!=helpState->state.end(); ++it ) {
  //   *it = proposal->Sample(state);
  // }

  // Run fused simulation
  problem->LogDensity(helpState);
  double* logDensityArray = boost::any_cast<double*>(helpState->meta["LogTarget"]);

  // Transfer LogDensity data to proposedStates
  unsigned int k = 0;
  for(auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) {
    (*it)->state[0] = helpState->state.at(k);
    (*it)->meta["LogTarget"] = logDensityArray[k++];
  }

  // evaluate the target density
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  for( unsigned int i=0; i<Np1; ++i )
    R(i) = AnyCast(proposedStates[i]->meta["LogTarget"]);

  // compute stationary transition probability
  AcceptanceDensity(R);
}