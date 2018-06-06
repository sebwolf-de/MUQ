#ifndef CRANKNICOLSONPROPOSAL_H_
#define CRANKNICOLSONPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    class CrankNicolsonProposal : public MCMCProposal {
    public:

      CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                            std::shared_ptr<AbstractSamplingProblem> prob,
                            std::shared_ptr<muq::Modeling::Gaussian> prior);

      CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                            std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~CrankNicolsonProposal() = default;

    protected:

      double beta;

      const Eigen::VectorXd priorMu;

      /// The proposal distribution
      std::shared_ptr<muq::Modeling::Gaussian> propPart;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) override;

      static std::shared_ptr<muq::Modeling::Gaussian> ExtractPrior(std::shared_ptr<AbstractSamplingProblem> prob, std::string nodeName);
    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
