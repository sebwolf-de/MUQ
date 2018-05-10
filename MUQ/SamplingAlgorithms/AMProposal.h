#ifndef AMPROPOSAL_H_
#define AMPROPOSAL_H_

#include "MUQ/SamplingAlgorithms/MHProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    class AMProposal : public MHProposal {
    public:

      AMProposal(boost::property_tree::ptree const& pt,
                 std::shared_ptr<AbstractSamplingProblem> prob);

      ~AMProposal();

      /// Adapt the proposal after each step
      /**
	 Adapt the proposal covariance.
	 @param[in] t The current step
	 @param[in] state The current state
       */
      virtual void Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states) override;

    private:

      /// Update the covariance of the samples
      /**
	 Adapt the proposal covariance.
	 @param[in] t The current step
	 @param[in] state The current state
       */
      void Update(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states);
      void UpdateOne(unsigned int const numSamps, std::shared_ptr<SamplingState> state);

      /// The current mean
      Eigen::VectorXd mean;

      /// The current covariance
      Eigen::MatrixXd cov;

      /// How frequently should we update the adaption?
      const unsigned int adaptSteps;

      /// When should we start adapting?
      const unsigned int adaptStart;

      // Multiplicative scaling of the sample covariance
      const double adaptScale;

    };

  } // namespace SamplingAlgorithms
} // namespace muq

#endif
