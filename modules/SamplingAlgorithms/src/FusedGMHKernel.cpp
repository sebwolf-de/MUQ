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
  std::cout << "Fused line 30" << std::endl;
  // If the current state does not have LogTarget information, add it
  if(! state->HasMeta("LogTarget"))
    state->meta["LogTarget"] = 0.0; // dummy value to avoid unfuded sim ... problem->LogDensity(state);

  std::shared_ptr<SamplingState> helpState = state;
  helpState->state.resize(N); // TODO: check for correct pointer syntax
  
  // propose the points
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;
  
  for(unsigned int j = 0; j<N; j++) {
    helpState->state.at(j) = proposal->Sample(state)->state[0];
  }
  // std::cout << "Size of N: " << N << std::endl;
  // std::cout << "Size of state vector: " << helpState->state.size() << std::endl;
  // for(auto it = helpState->state.begin()+1; it!=helpState->state.end(); ++it ) {
  //   *it = proposal->Sample(state);
  // }
  std::cout << "Fused line 50" << std::endl;
  // Run fused simulation
  problem->LogDensity(helpState);
  double* logDensityArray = boost::any_cast<double*>(helpState->meta["LogTarget"]);
  std::cout << "Fused line 54" << std::endl;
  // Transfer LogDensity data to proposedStates
  unsigned int k = 0;
  for(auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) {
    *it = helpState->state.at(k);
    (*it)->meta["LogTarget"] = logDensityArray[k++];
  }
  std::cout << "Fused line 61" << std::endl;
  // evaluate the target density
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  for( unsigned int i=0; i<Np1; ++i )
    R(i) = AnyCast(proposedStates[i]->meta["LogTarget"]);

  std::cout << "Fused line 67" << std::endl;
  // compute stationary transition probability
  AcceptanceDensity(R);
}