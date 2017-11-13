#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

//#inclue <mpi.h>

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace SamplingAlgorithms {
    
    class SamplingAlgorithm : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm();
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
