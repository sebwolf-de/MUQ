#ifndef MHKERNEL_H_
#define MHKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace SamplingAlgorithms {
    /// Monte Carlo transition kernel
    /**
       Samples from the target distirbution directly and returns that state.
     */
    class MHKernel : public TransitionKernel {
    public:

      MHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      ~MHKernel();

      virtual std::shared_ptr<SamplingState> Step(std::shared_ptr<SamplingState> prevState) override;

      // What block of the state does this kernel work on?
      const int blockInd;

    protected:
      std::shared_ptr<muq::Modeling::Distribution> proposal;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
