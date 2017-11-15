#ifndef MHKERNEL_H_
#define MHKERNEL_H_

#include "MUQ/SamplingAlgorithms/MCMCKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
    /// Monte Carlo transition kernel
    /**
       Samples from the target distirbution directly and returns that state.
     */
    class MHKernel : public MCMCKernel {
    public:

      MHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      ~MHKernel();
      
    private:
      
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
