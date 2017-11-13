#ifndef SAMPLINGPROBLEM_H_
#define SAMPLINGPROBLEM_H_

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace SamplingAlgorithms {
    class SamplingProblem {
    public:

      SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> target);

      ~SamplingProblem();

      /// Directly sample the target distribution
      /**
	 Assumes that muq::SamplingAlgorithms::SamplingProblem::target has a Sample method implemented.
	 @param[in] inputs Inputs to the target distribution Sample method
	 \return A sample from the target distribution
       */
      boost::any SampleTarget(muq::Modeling::ref_vector<boost::any> const& inputs) const;
      
    private:

      /// The target distribution (the prior in the inference case)
      std::shared_ptr<muq::Modeling::Distribution> target;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
