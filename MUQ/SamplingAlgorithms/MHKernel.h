#ifndef MHKERNEL_H_
#define MHKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

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
      
    private:
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
