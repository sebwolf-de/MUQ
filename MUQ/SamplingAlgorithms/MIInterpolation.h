#ifndef MIINTERPOLATION_H_
#define MIINTERPOLATION_H_

namespace muq {
  namespace SamplingAlgorithms {

    class MIInterpolation {
    public:

      //MIInterpolation() {}

      //~MIInterpolation() {}

      virtual std::shared_ptr<SamplingState> interpolate (std::shared_ptr<SamplingState> coarseProposal, std::shared_ptr<SamplingState> fineProposal) = 0;


    };


  }
}

#endif
