#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

MCMCProposal::MCMCProposal() : Distribution() {}

MCMCProposal::~MCMCProposal() {}

std::shared_ptr<MCMCProposal> MCMCProposal::Construct(pt::ptree const& pt) {
    // get the name of the proposal
  const std::string& proposalName = pt.get<std::string>("MCMC.Proposal");

  // construct it from the map
  return GetMCMCProposalMap()->at(proposalName) (pt);  
}

std::shared_ptr<MCMCProposal::MCMCProposalMap> MCMCProposal::GetMCMCProposalMap() {
  // define a static map from type to constructor
  static std::shared_ptr<MCMCProposalMap> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<MCMCProposalMap>();
  }

  return map;
}
