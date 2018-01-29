#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_MCMC_PROPOSAL(AMProposal)
AMProposal::AMProposal(pt::ptree const& pt) : MHProposal(pt),
                                              adaptSteps(pt.get<unsigned int>("MCMC.Proposal.AM.Steps", 1)),
                                              adaptStart(pt.get<unsigned int>("MCMC.Proposal.AM.Start", 1)) {}


AMProposal::~AMProposal() {}

void AMProposal::Adapt(unsigned int const t, std::shared_ptr<SamplingState> state) {
  // always update the sample mean and covariance
  Update(t, state);

  if( t%adaptSteps==0 && t>adaptStart ) {
    // the new proposal covariance
    boost::any propCov = algebra->Identity(cov.type(), algebra->Size(cov, 0), algebra->Size(cov, 1));
    propCov = algebra->Add(algebra->Multiply(1.0, cov), algebra->Multiply(1.0e-10, propCov));

    // update the proposal covariance
    proposal->SetCovariance(propCov);
  }
}

void AMProposal::Update(unsigned int const t, std::shared_ptr<SamplingState> state) {
  // first sample---we have no mean, just set it to the first sample
  if( mean.type()==typeid(boost::none) ) { mean = state->state; return; }

  // need to store the old mean
  const boost::any oldMean = mean;

  // update the mean
  const double dum = 1.0/(1.0+(double)t);
  mean = algebra->Add(algebra->Multiply((double)t*dum, mean), algebra->Multiply(dum, state->state));

  // second sample---we have a mean, but not cov
  if( cov.type()==typeid(boost::none) ) { cov = algebra->OuterProduct(algebra->Subtract(state->state, oldMean), algebra->Subtract(state->state, mean)); return; }

  // update the cov
  cov = algebra->Add(algebra->Multiply((double)t*dum, cov), algebra->Multiply(1.0/(double)t, algebra->OuterProduct(algebra->Subtract(state->state, oldMean), algebra->Subtract(state->state, mean))));
}
