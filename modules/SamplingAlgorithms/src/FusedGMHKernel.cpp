#include "MUQ/SamplingAlgorithms/FusedGMHKernel.h"

#include <Eigen/Eigenvalues>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/AnyHelpers.h"

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
  // If the current state does not have LogTarget information, evaluate and add it
  // Should only happen once in the beginning.
  if(! state->HasMeta("LogTarget")) {
    std::cout << "State does not have meta info, run Sampling Problem." << std::endl;
    std::shared_ptr<SamplingState> helpState = std::make_shared<SamplingState>(state->state);
    helpState->state.resize(N);
    for (size_t i = 0; i < N; i++) {
      helpState->state.at(i) = state->state[0];
    }
    problem->LogDensity(helpState);
    std::vector<double> logDensityArray = boost::any_cast<std::vector<double>>(helpState->meta["LogTarget"]);
    state->meta["LogTarget"] = logDensityArray.at(0);
  }
    
  std::shared_ptr<SamplingState> helpState = std::make_shared<SamplingState>(state->state);

  helpState->state.resize(N);
  
  // propose the points
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;
  
  for(unsigned int j = 0; j<N; j++) {
    helpState->state.at(j) = proposal->Sample(state)->state[0];
  }
  
  // Run fused simulation
  problem->LogDensity(helpState);
  std::vector<double> logDensityArray = boost::any_cast<std::vector<double>>(helpState->meta["LogTarget"]);
  
  // Transfer LogDensity data to proposedStates
  unsigned int k = 0;
  for(auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) {
    *it = std::make_shared<SamplingState>(helpState->state.at(k));
    (*it)->meta["LogTarget"] = logDensityArray[k++];
  }
  
  // evaluate the target density
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  R(0) = boost::any_cast<double>(state->meta["LogTarget"]);
  for( unsigned int i=1; i<Np1; ++i ) {
    R(i) = logDensityArray[i-1];
  }
  
  // compute stationary transition probability
  AcceptanceDensity(R);
}
