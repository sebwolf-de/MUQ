#ifndef MHPROPOSAL_H_
#define MHPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    class MHProposal : public MCMCProposal {
    public:

      MHProposal(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~MHProposal() = default;

    protected:

      /// The proposal distribution
      std::shared_ptr<muq::Modeling::Gaussian> proposal;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) override;


    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
