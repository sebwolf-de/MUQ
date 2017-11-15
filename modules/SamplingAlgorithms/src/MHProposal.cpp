#include "MUQ/SamplingAlgorithms/MHProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_MCMC_PROPOSAL(MHProposal)

MHProposal::MHProposal(pt::ptree const& pt) : MCMCProposal() {
  // compute the (diagonal) covariance for the proposal
  const Eigen::VectorXd cov = pt.get<double>("MCMC.Proposal.MH.ProposalSize", 1.0)*Eigen::VectorXd::Ones(pt.get<unsigned int>("MCMC.StateDimension"));
  
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
  boost::any prop = proposal->Sample(std::pair<boost::any, Gaussian::Mode>(current->state, Gaussian::Mode::Mean));

  // store the new state in the output
  return std::make_shared<SamplingState>(prop, 1.0);
}

double MHProposal::LogDensityImpl(ref_vector<boost::any> const& inputs) {
  // get the state
  std::shared_ptr<SamplingState> state = boost::any_cast<std::shared_ptr<SamplingState> >(inputs[0]);
  assert(state);

  // get the conditioned state
  std::shared_ptr<SamplingState> conditioned = boost::any_cast<std::shared_ptr<SamplingState> >(inputs[1]);
  assert(conditioned);

  return proposal->LogDensity(state->state, std::pair<boost::any, Gaussian::Mode>(conditioned->state, Gaussian::Mode::Mean));
}
