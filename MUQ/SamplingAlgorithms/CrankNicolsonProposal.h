#ifndef CRANKNICOLSONPROPOSAL_H_
#define CRANKNICOLSONPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/ModPiece.h"

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

      // When the prior distribution has inputs, we need to save information to evaluate them
      std::shared_ptr<muq::Modeling::ModPiece> priorMeanModel;
      std::vector<int> priorMeanInds;

      std::shared_ptr<muq::Modeling::ModPiece> priorCovModel;
      std::vector<int> priorCovInds;
      bool priorUsesCov;

      std::vector<Eigen::VectorXd>GetPriorInputs(std::vector<Eigen::VectorXd> const& currState);

      //const Eigen::VectorXd priorMu;

      /// The proposal distribution
      //std::shared_ptr<muq::Modeling::Gaussian> propPart;

      // Sometimes we have to keep track of the prior distribution so we can update the proposal mean and covariance
      std::shared_ptr<muq::Modeling::Gaussian> priorDist;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) override;
      
      void ExtractPrior(std::shared_ptr<AbstractSamplingProblem> prob, std::string nodeName);
    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
