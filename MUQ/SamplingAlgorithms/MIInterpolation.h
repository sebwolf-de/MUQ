#ifndef MIINTERPOLATION_H_
#define MIINTERPOLATION_H_

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace muq {
  namespace SamplingAlgorithms {

    class MIInterpolation {
    public:

      virtual ~MIInterpolation() = default;

      virtual std::shared_ptr<SamplingState> interpolate (std::shared_ptr<SamplingState> coarseProposal, std::shared_ptr<SamplingState> fineProposal) = 0;


    };


  }
}

#endif
