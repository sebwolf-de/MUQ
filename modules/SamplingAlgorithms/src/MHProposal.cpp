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

boost::any MHProposal::SampleImpl(ref_vector<boost::any> const& inputs) {
  // get the current state
  std::shared_ptr<SamplingState> current = boost::any_cast<std::shared_ptr<SamplingState> >(inputs[0]);
  assert(current);

  // the mean of the proposal is the current point
  std::vector<boost::any> props = current->state;
  Eigen::VectorXd const& xc = AnyConstCast(current->state.at(blockInd));

  boost::any anyProp = proposal->AsVariable()->Sample();
  Eigen::VectorXd const& prop = AnyConstCast(anyProp);

  props.at(blockInd) =  (xc + prop).eval();

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double MHProposal::LogDensityImpl(ref_vector<boost::any> const& inputs) {
  // get the state
  std::shared_ptr<SamplingState> state = boost::any_cast<std::shared_ptr<SamplingState> >(inputs[0]);
  assert(state);

  // get the conditioned state
  std::shared_ptr<SamplingState> conditioned = boost::any_cast<std::shared_ptr<SamplingState> >(inputs[1]);
  assert(conditioned);

  Eigen::VectorXd const& a = AnyConstCast(state->state.at(blockInd));
  Eigen::VectorXd const& b = AnyConstCast(conditioned->state.at(blockInd));

  return proposal->LogDensity(boost::any((b-a).eval()));//, std::pair<boost::any, Gaussian::Mode>(conditioned->state.at(blockInd), Gaussian::Mode::Mean));
}
