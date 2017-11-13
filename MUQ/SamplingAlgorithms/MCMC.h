#ifndef MCMC_H_
#define MCMC_H_

#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace muq {
  namespace SamplingAlgorithms {
    class MCMC : public SamplingAlgorithm {
    public:

      MCMC();
      
      ~MCMC();
      
    private:

      /// Create the transition kernel
      /**
	 @param[in] pt Parameters for the kernel
	 @param[in] problem The sampling problem that computes the next state in the MCMC chain
	 \return The transition kernel
       */
      virtual std::shared_ptr<TransitionKernel> Kernel(boost::property_tree::ptree& pt, std::shared_ptr<SamplingProblem> problem) const override;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
