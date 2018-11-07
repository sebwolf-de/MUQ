#ifndef MIINTERPOLATION_H_
#define MIINTERPOLATION_H_

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Interpolation interface combining coarse and fine samples.
        @details This interface defines how coarse and fine samples are
        combined, as needed for Multiindex MCMC proposals.
     */
    class MIInterpolation {
    public:

      virtual ~MIInterpolation() = default;

      virtual std::shared_ptr<SamplingState> interpolate (std::shared_ptr<SamplingState> coarseProposal, std::shared_ptr<SamplingState> fineProposal) = 0;


    };


  }
}

#endif
