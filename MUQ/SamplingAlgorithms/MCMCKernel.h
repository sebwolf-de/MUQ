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

      /// Allow the kernel to adapt given a new state
      /**
	 By default this function does nothing but children can override it to adapt the kernel
	 @param[in] t The current step
	 @param[in] state The current state
       */
      virtual void PostStep(unsigned int const t, std::shared_ptr<SamplingState> state) override;
      
    protected:

      /// The proposal
      std::shared_ptr<MCMCProposal> proposal;
      
    private:
      
    };
    
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
