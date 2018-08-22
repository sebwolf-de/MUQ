#ifndef PARALLELAMPROPOSAL_H_
#define PARALLELAMPROPOSAL_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER

#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ParallelAMProposal : public AMProposal {
    public:

      ParallelAMProposal(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem);

      ~ParallelAMProposal() = default;

    private:
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // end MUQ_HAS_PARCER
#endif
