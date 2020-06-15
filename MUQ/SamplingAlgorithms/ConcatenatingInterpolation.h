#ifndef CONCATENATINGINTERPOLATION_H_
#define CONCATENATINGINTERPOLATION_H_

#include "MUQ/SamplingAlgorithms/MIInterpolation.h"
#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

namespace muq {
  namespace SamplingAlgorithms {

    class ConcatenatingInterpolation : public MIInterpolation {
    public:
      ConcatenatingInterpolation(std::shared_ptr<muq::Utilities::MultiIndex> const& index) : index(index) {
    	}

    	virtual std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) override {
    		int fine_part_size = fineProposal->state[0].size() - coarseProposal->state[0].size();

    		Eigen::VectorXd interpolatedState(fineProposal->state[0].size());
    		interpolatedState << coarseProposal->state[0], fineProposal->state[0].tail(fine_part_size);

    		return std::make_shared<SamplingState>(interpolatedState);
    	}

    private:
      std::shared_ptr<muq::Utilities::MultiIndex> index;
    };

  }
}


#endif
