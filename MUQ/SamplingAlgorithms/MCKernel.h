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

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// The number of Monte Carlo samples
      const unsigned int N;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
