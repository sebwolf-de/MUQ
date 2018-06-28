#ifndef SAMPLINGPROBLEM_H_
#define SAMPLINGPROBLEM_H_

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
    @ingroup SamplingAlgorithms
    @class SamplingProblem
    @brief Class for sampling problems based purely on a density function.
    */
    class SamplingProblem : public AbstractSamplingProblem{
    public:

      /**
	     @param[in] target The target distribution
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> targetIn);

      virtual ~SamplingProblem() = default;


      virtual double LogDensity(std::shared_ptr<SamplingState> state) override;

      virtual Eigen::VectorXd GradLogDensity(std::shared_ptr<SamplingState> state,
                                             unsigned                       blockWrt) override;


      std::shared_ptr<muq::Modeling::ModPiece> GetDistribution(){return target;};

    private:

      /// The target distribution (the prior in the inference case)
      std::shared_ptr<muq::Modeling::ModPiece> target;

      static unsigned GetNumBlocks(std::shared_ptr<muq::Modeling::ModPiece> target);
      static std::vector<int> GetBlockSizes(std::shared_ptr<muq::Modeling::ModPiece> target);

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
