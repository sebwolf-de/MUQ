#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace SamplingAlgorithms {
    
    class SamplingAlgorithm : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm();
      
    private:

      /**
	 Inputs:
	 <ol>
	 <li> Parameters for the algorithm
	 </ol>
	 @param[in] inputs Inputs to the algorithm
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
