#ifndef MCKERNEL_H_
#define MCKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
    /// Monte Carlo transition kernel
    /**
       Samples from the target distirbution directly and returns that state.
     */
    class MCKernel : public TransitionKernel {
    public:

      MCKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      ~MCKernel();
      
    private:
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
