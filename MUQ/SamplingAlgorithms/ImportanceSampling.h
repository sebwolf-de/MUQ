#ifndef IMPORTANCESAMPLING_H_
#define IMPORTANCESAMPLING_H_

#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ImportanceSampling : public SamplingAlgorithm {
    public:

      ImportanceSampling();

      ~ImportanceSampling();
      
    private:

      /// Create the transition kernel
      /**
	 @param[in] pt Parameters for the kernel
	 @param[in] problem The sampling problem that evaluates the distribution we are trying to characterize and samples the biasing distribution
	 \return The transition kernel
       */
      virtual std::shared_ptr<TransitionKernel> Kernel(boost::property_tree::ptree& pt, std::shared_ptr<SamplingProblem> problem) const override;

    };
  } // namespace SamplingAlgorithms
} // namespace muq


#endif
