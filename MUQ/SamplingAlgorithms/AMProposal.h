#ifndef AMPROPOSAL_H_
#define AMPROPOSAL_H_

#include "MUQ/SamplingAlgorithms/MHProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    class AMProposal : public MHProposal {
    public:

      AMProposal(boost::property_tree::ptree const& pt);

      ~AMProposal();

      /// Adapt the proposal after each step
      /**
	 Adapt the proposal covariance.
	 @param[in] t The current step
	 @param[in] state The current state
       */
      virtual void Adapt(unsigned int const t, std::shared_ptr<SamplingState> state) override;
      
    private:

      /// Update the covariance of the samples
      /**
	 Adapt the proposal covariance.
	 @param[in] t The current step
	 @param[in] state The current state
       */
      void Update(unsigned int const t, std::shared_ptr<SamplingState> state);

      /// The current mean
      boost::any mean = boost::none;

      /// The current covariance
      boost::any cov = boost::none;

      /// How frequently should we update the adaption?
      const unsigned int adaptSteps;

      /// When should we start adapting?
      const unsigned int adaptStart;

    };
    
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
