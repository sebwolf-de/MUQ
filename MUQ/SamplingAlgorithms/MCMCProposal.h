#ifndef MCMCPROPOSAL_H_
#define MCMCPROPOSAL_H_

#include <map>

#include <boost/function.hpp>
#include <boost/property_tree/ptree.hpp>

#include "MUQ/Utilities/RegisterClassName.h"

#include "MUQ/Modeling/Distributions/Distribution.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace muq {
  namespace SamplingAlgorithms {

    class MCMCProposal : public muq::Modeling::Distribution {
    public:

      MCMCProposal();

      ~MCMCProposal();

      /// Static constructor for the transition kernel
      /**
	 @param[in] pt The options for the MCMC kernel
	 \return The MCMC proposal
       */
      static std::shared_ptr<MCMCProposal> Construct(boost::property_tree::ptree const& pt);

      typedef boost::function<std::shared_ptr<MCMCProposal>(boost::property_tree::ptree)> MCMCProposalConstructor;

      typedef std::map<std::string, MCMCProposalConstructor> MCMCProposalMap;

      static std::shared_ptr<MCMCProposalMap> GetMCMCProposalMap();
      
    private:

    };    
  } // namespace SamplingAlgoirthms
} // namespace muq

#define REGISTER_MCMC_PROPOSAL(NAME) static auto reg ##NAME		\
  = muq::SamplingAlgorithms::MCMCProposal::GetMCMCProposalMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));
    

#endif
