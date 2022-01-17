#ifndef FusedGMHKERNEL_H_
#define FusedGMHKERNEL_H_

#include "MUQ/config.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/GMHKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
      /// A kernel for the generalized Metropolis-Hastings kernel
    /**
       Reference: "A general construction for parallelizing Metropolis-Hastings algorithms" (Calderhead, 2014)
     */
    class FusedGMHKernel : public GMHKernel {
    public:
    /**
	 @param[in] pt Options for this kenel and the standard Metropolis-Hastings kernel
	 @param[in] problem The problem we want to sample
       */
      FusedGMHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem);

      /**
	 @param[in] pt Options for this kenel and the standard Metropolis-Hastings kernel
	 @param[in] problem The problem we want to sample
	 @param[in] proposalIn The proposal for the MCMC chain
       */
      FusedGMHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn);

      virtual ~FusedGMHKernel();

    private:

    /**
	 Propose GMHKernel::N points and compute the cumulative distribution of the stationary distribution for the acceptance probability (GMHKernel::proposedStates and GMHKernel::stationaryAcceptance, respectively)
	 @param[in] t The current step in the MCMC chain
	 @param[in] state The current MCMC state
      */
    virtual void PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) override;

    /// Adapt the function to a single fused simulation call instead of multiple serial ones
    /// Propose \f$N\f$ points in serial and evaluate the log target
    /**
    @param[in] t The current step in the MCMC chain
    @param[in] state The current point
    */
    void FusedProposal(unsigned int const t, std::shared_ptr<SamplingState> state);


    };



  } // namespace SamplingAlgorithms
} // namespace muq

#endif