#ifndef MHKERNEL_H_
#define MHKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {
    /// Monte Carlo transition kernel
    /**
       Samples from the target distirbution directly and returns that state.
     */
    class MHKernel : public TransitionKernel {
    public:

      MHKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem);

      MHKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem,
               std::shared_ptr<MCMCProposal> proposalIn);

      virtual ~MHKernel() = default;

      virtual inline std::shared_ptr<MCMCProposal> Proposal() {return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

      virtual void PrintStatus(std::string prefix) const override;


      virtual inline double AcceptanceRate() const {return double(numAccepts)/double(numCalls);};

#if MUQ_HAS_PARCER
      virtual void SetCommunicator(std::shared_ptr<parcer::Communicator> newcomm) override;
#endif

    protected:
      std::shared_ptr<MCMCProposal> proposal;

      unsigned int numCalls = 0;
      unsigned int numAccepts = 0;

      /// true: reevaluate the log density (even if one already exists), false: use stored log density
      /**
	 For example, if the log-density is a continually refined surrogate (LA-MCMC) or an importance sampling estimate (pseudo-marginal) then we need to reevaluate every time.
       */
      const bool reeval;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
