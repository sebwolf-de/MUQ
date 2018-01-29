#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

MCMCProposal::MCMCProposal() : Distribution() {}

MCMCProposal::~MCMCProposal() {}

std::shared_ptr<MCMCProposal> MCMCProposal::Construct(pt::ptree const& pt) {

  // get the name of the proposal
  std::string proposalName = pt.get("MCMC.Proposal", "MHProposal");

  auto proposalMap = GetMCMCProposalMap();
  auto iter = proposalMap->find(proposalName);

  if(iter == proposalMap->end()){
    std::cerr << "ERROR: Could not find MCMC proposal \"" << proposalName << "\".  Available options are:\n";

    for(auto it=proposalMap->begin(); it!=proposalMap->end(); ++it)
      std::cerr << "  " << it->first << std::endl;
    std::cerr << std::endl;

    assert(iter!=proposalMap->end());
  }

  return iter->second(pt);

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

void MCMCProposal::Adapt(unsigned int const t, std::shared_ptr<SamplingState> state) {}
