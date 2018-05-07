#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/Utilities/AnyHelpers.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(AMProposal)
AMProposal::AMProposal(pt::ptree const& pt , std::shared_ptr<AbstractSamplingProblem> prob) : MHProposal(pt, prob),
                                              adaptSteps(pt.get<unsigned int>("AdaptSteps")),
                                              adaptStart(pt.get<unsigned int>("AdaptStart")),
                                              adaptScale(pt.get<double>("AdaptScale")) {}


AMProposal::~AMProposal() {}

void AMProposal::Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states) {

  // always update the sample mean and covariance
  Update(t, states);

  if( t%adaptSteps==0 && t>adaptStart ) {
    // the new proposal covariance
    Eigen::MatrixXd const& sampCov = AnyConstCast(cov);

    Eigen::MatrixXd adjustedCov = adaptScale * sampCov + 1e-10 * Eigen::MatrixXd::Identity(sampCov.rows(), sampCov.cols());

    // update the proposal covariance
    proposal->SetCovariance(adjustedCov);
  }
}

void AMProposal::UpdateOne(unsigned int const numSamps, std::shared_ptr<SamplingState> state)
{
  // first sample---we have no mean, just set it to the first sample
  if( mean.type()==typeid(boost::none) ){
    mean = state->state.at(blockInd);
    return;
  }

  // update the mean
  Eigen::VectorXd oldMean = AnyCast(mean);
  Eigen::VectorXd& newMean = AnyCast(mean);
  Eigen::VectorXd const& newState  = AnyConstCast(state);

  newMean = (oldMean*numSamps + newState)/(numSamps+1.0);

  // If we haven't compute the covariance before...
  if( cov.type()==typeid(boost::none) ){

    cov = Eigen::MatrixXd::Zero(oldMean.size(),oldMean.size());
    Eigen::MatrixXd& newCov = AnyCast(cov);

    //compute covariance from scratch, from the definition
    newCov.selfadjointView<Eigen::Lower>().rankUpdate(oldMean - newMean, 1.0);
    newCov.selfadjointView<Eigen::Lower>().rankUpdate(newState - newMean, 1.0);

  }else{
    Eigen::MatrixXd& newCov = AnyCast(cov);
    newCov *= (numSamps - 1.0) / numSamps;

    //note that the asymmetric form fixes the fact that the old mean was wrong
    newCov += (1.0 / static_cast<double>(numSamps)) * (newState - oldMean) * (newState - newMean).transpose();
  }
}

void AMProposal::Update(unsigned int const numSamps, std::vector<std::shared_ptr<SamplingState>> const& states) {

  for(int i=0; i<states.size(); ++i)
    UpdateOne(numSamps+i,states.at(i));

}
