#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"

#include "MUQ/Utilities/AnyHelpers.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(MHProposal)

MHProposal::MHProposal(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob) : MCMCProposal(pt,prob) {

  unsigned int problemDim = prob->blockSizes.at(0);

  // compute the (diagonal) covariance for the proposal
  const Eigen::VectorXd cov = pt.get("ProposalVariance", 1.0)*Eigen::VectorXd::Ones(problemDim);

  // created a Gaussian with scaled identity covariance
  if( cov.size()==1 ) {
    proposal = std::make_shared<muq::Modeling::Gaussian>(cov(0), Gaussian::Mode::Covariance);
  } else {
    proposal = std::make_shared<muq::Modeling::Gaussian>(cov, Gaussian::Mode::Covariance);
  }
}

MHProposal::~MHProposal() {}

std::shared_ptr<SamplingState> MHProposal::Sample(std::shared_ptr<SamplingState> currentState) {

  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;
  Eigen::VectorXd const& xc = currentState->state.at(blockInd);

  boost::any anyProp = proposal->AsVariable()->Sample();
  Eigen::VectorXd const& prop = AnyConstCast(anyProp);

  props.at(blockInd) = xc + prop;

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double MHProposal::LogDensity(std::shared_ptr<SamplingState> currState,
                              std::shared_ptr<SamplingState> propState) {

  Eigen::VectorXd diff = currState->state.at(blockInd)-propState->state.at(blockInd);
  return proposal->LogDensity(boost::any(diff));//, std::pair<boost::any, Gaussian::Mode>(conditioned->state.at(blockInd), Gaussian::Mode::Mean));
}
