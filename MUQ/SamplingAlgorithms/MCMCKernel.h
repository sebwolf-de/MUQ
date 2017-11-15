#ifndef MCMCKERNEL_H_
#define MCMCKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    class MCMCKernel : public TransitionKernel {
    public:

      MCMCKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      ~MCMCKernel();
      
    protected:

      /// The proposal
      std::shared_ptr<MCMCProposal> proposal;
      
    private:
      
    };
    
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
