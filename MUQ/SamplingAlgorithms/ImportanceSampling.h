#ifndef IMPORTANCESAMPLING_H_
#define IMPORTANCESAMPLING_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Distribution.h"

#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace muq {
  namespace SamplingAlgorithms {
    /** @ingroup MCMC
      @class ImportanceSampling
      @brief Defines an imporance sampline sampler
      @details
      <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "NumSamples"  | Int | - | The total number of steps (including burnin) to take, i.e., the length of the Markov chain. |
    */
    class ImportanceSampling : public SamplingAlgorithm {
    public:

      ImportanceSampling(std::shared_ptr<muq::Modeling::Density> const& target, std::shared_ptr<muq::Modeling::Distribution> const& bias, boost::property_tree::ptree const& pt);

      ~ImportanceSampling();

    private:

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

      /// The number of samples
      const unsigned int numSamps;

      // The target density
      std::shared_ptr<muq::Modeling::Density> target;

      /// The biasing distribution
      std::shared_ptr<muq::Modeling::Distribution> bias;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
