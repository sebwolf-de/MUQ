#ifndef SAMPLINGPROBLEM_H_
#define SAMPLINGPROBLEM_H_

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace SamplingAlgorithms {
    class SamplingProblem {
    public:

      /**
	 @param[in] target The target distribution
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> target);

      /**
	 @param[in] target The target distribution
	 @param[in] bias A biasing distribution
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> target, std::shared_ptr<muq::Modeling::Distribution> bias);

      ~SamplingProblem();

      /// Directly sample the target distribution
      /**
	 Assumes that muq::SamplingAlgorithms::SamplingProblem::target has a Sample method implemented.
	 @param[in] inputs Inputs to the target distribution Sample method
	 \return A sample from the target distribution
       */
      boost::any SampleTarget(muq::Modeling::ref_vector<boost::any> const& inputs) const;

      /// Evaluate the log target distribution
      /**
	 Assumes that muq::SamplingAlgorithms::SamplingProblem::target has a LogDensity method implemented.
	 @param[in] inputs Inputs to the target distribution LogDensity method
	 \return The log density
       */
      double EvaluateLogTarget(muq::Modeling::ref_vector<boost::any> const& inputs) const;

      /// Sample the biasing distribution
      /**
	 Assumes that muq::SamplingAlgorithms::SamplingProblem::bias has a Sample method implemented.
	 @param[in] inputs Inputs to the biasing distribution Sample method
	 \return A sample from the biasing distribution
       */
      boost::any SampleBiasingDistribution(muq::Modeling::ref_vector<boost::any> const& inputs) const;

      /// Evaluate the log biasing distribution
      /**
	 Assumes that muq::SamplingAlgorithms::SamplingProblem::bias has a LogDensity method implemented.
	 @param[in] inputs Inputs to the biasing distribution LogDensity method
	 \return The log density
       */
      double EvaluateLogBiasingDistribution(muq::Modeling::ref_vector<boost::any> const& inputs) const;

    private:

      /// The target distribution (the prior in the inference case)
      std::shared_ptr<muq::Modeling::Distribution> target;

      /// An optional biasing distribution
      boost::optional<std::shared_ptr<muq::Modeling::Distribution> > bias = boost::none;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
