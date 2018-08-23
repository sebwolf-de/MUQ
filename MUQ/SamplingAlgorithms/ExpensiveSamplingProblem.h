#ifndef EXPENSIVESAMPLINGPROBLEM_H_
#define EXPENSIVESAMPLINGPROBLEM_H_

#include "MUQ/Approximation/Regression/LocalRegression.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ExpensiveSamplingProblem : public SamplingProblem {
    public:

      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree const& pt);

      ~ExpensiveSamplingProblem() = default;

      virtual double LogDensity(std::shared_ptr<SamplingState> state) override;
      
    private:

      std::shared_ptr<muq::Approximation::LocalRegression> reg;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
