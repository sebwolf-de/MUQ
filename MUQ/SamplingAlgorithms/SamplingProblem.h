#ifndef SAMPLINGPROBLEM_H_
#define SAMPLINGPROBLEM_H_

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Class for sampling problems based purely on a density function.
    */
    class SamplingProblem : public AbstractSamplingProblem{
    public:

      SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> targetIn,
                      std::vector<int>                      const& inputSizes);
      /**
	     @param[in] target The target distribution
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> targetIn);

      virtual ~SamplingProblem() = default;


      virtual double LogDensity(std::shared_ptr<SamplingState> state) override;

      virtual boost::any GradLogDensity(std::shared_ptr<SamplingState> state,
                                        unsigned                       blockWrt) override;


      std::shared_ptr<muq::Modeling::Distribution> GetDistribution(){return target;};

    private:

      /// The target distribution (the prior in the inference case)
      std::shared_ptr<muq::Modeling::Distribution> target;

      static unsigned GetNumBlocks(std::shared_ptr<muq::Modeling::Distribution> target);
      static std::vector<int> GetBlockSizes(std::shared_ptr<muq::Modeling::Distribution> target);

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
