#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Modeling;

REGISTER_MCMC_PROPOSAL(CrankNicolsonProposal)

CrankNicolsonProposal::CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                                             std::shared_ptr<AbstractSamplingProblem> prob,
                                             std::shared_ptr<Gaussian>                prior) : MCMCProposal(pt,prob),
                                                                                               beta(pt.get("Beta",0.5)),
                                                                                               priorMu(prior->GetMean())
{
  int dim = priorMu.size();

  // Construct the random part of the proposal from the prior covariance
  if(prior->GetMode() == Gaussian::Covariance){
    propPart = std::make_shared<Gaussian>(Eigen::VectorXd::Zero(dim), prior->GetCovariance(), Gaussian::Covariance );
  }else{
    propPart = std::make_shared<Gaussian>(Eigen::VectorXd::Zero(dim), prior->GetPrecision(), Gaussian::Precision );
  }
}

CrankNicolsonProposal::CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                                             std::shared_ptr<AbstractSamplingProblem> prob) : CrankNicolsonProposal(pt,prob,ExtractPrior(prob, pt.get<std::string>("PriorNode")))
{}


std::shared_ptr<SamplingState> CrankNicolsonProposal::Sample(std::shared_ptr<SamplingState> currentState)
{
  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;

  props.at(blockInd) = priorMu + sqrt(1.0-beta*beta)*(currentState->state.at(blockInd)-priorMu) + beta*propPart->Sample();

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double CrankNicolsonProposal::LogDensity(std::shared_ptr<SamplingState> currState,
                                         std::shared_ptr<SamplingState> propState)
{
  Eigen::VectorXd diff = (propState->state.at(blockInd) - priorMu - sqrt(1.0-beta*beta)*(currState->state.at(blockInd)-priorMu))/beta;

  return propPart->LogDensity(diff);
}

std::shared_ptr<Gaussian> CrankNicolsonProposal::ExtractPrior(std::shared_ptr<AbstractSamplingProblem> prob,
                                                              std::string                              nodeName)
{
  // Cast the abstract base class into a sampling problem
  std::shared_ptr<SamplingProblem> prob2 = std::dynamic_pointer_cast<SamplingProblem>(prob);
  assert(prob2);

  // From the sampling problem, extract the ModPiece and try to cast it to a ModGraphPiece
  std::shared_ptr<ModPiece> targetDens = prob2->GetDistribution();
  std::shared_ptr<ModGraphPiece> targetDens2 = std::dynamic_pointer_cast<ModGraphPiece>(targetDens);
  assert(targetDens2);

  // Get the graph
  auto graph = targetDens2->GetGraph();

  // Get the prior piece corresponding to the Gaussian name
  auto priorPiece = graph->GetPiece(nodeName);
  assert(priorPiece);

  // Get the prior distribution
  std::shared_ptr<Density> priorDens = std::dynamic_pointer_cast<Density>(priorPiece);
  assert(priorDens);

  std::shared_ptr<Gaussian> gaussDist = std::dynamic_pointer_cast<Gaussian>(priorDens->GetDistribution());
  assert(gaussDist);

  return gaussDist;
}
