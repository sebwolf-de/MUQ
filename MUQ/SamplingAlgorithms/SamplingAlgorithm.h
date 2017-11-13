#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
    
    class SamplingAlgorithm : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm();

      ~SamplingAlgorithm();
      
    private:

      /**
	 Inputs:
	 <ol>
	 <li> Parameters for the algorithm
	 <li> The muq::SamplingAlgorithms::SamplingProblem that evaluates/samples the target distribution
	 </ol>
	 @param[in] inputs Inputs to the algorithm
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// Generate the next sample
      /**
	 @param[in] problem The muq::SamplingAlgorithms::SamplingProblem that evaluates/samples the target distribution
       */
      //virtual boost::any SampleOnce(std::shared_ptr<SamplingProblem> problem) const;

      /// Create the transition kernel
      /**
	 @param[in] pt Parameters for the kernel
	 @param[in] problem The sampling problem that evaluates/samples the distribution we are trying to characterize
	 \return The transition kernel
       */
      virtual std::shared_ptr<TransitionKernel> Kernel(boost::property_tree::ptree& pt, std::shared_ptr<SamplingProblem> problem) const = 0;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
