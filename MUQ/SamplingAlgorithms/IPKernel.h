#ifndef IPKERNEL_H_
#define IPKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

namespace muq {
  namespace SamplingAlgorithms {

    /// Importance sampling transition kernel
    /**
       Propose from a biasing distribution, compute the weight and return a state with those values
     */
    class IPKernel : public TransitionKernel {
    public:

      IPKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      ~IPKernel();

    private:
    };
    
  } // namespace SamplingAlgorithms
} // namespace muq


#endif
