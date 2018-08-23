#ifndef EXPENSIVESAMPLINGPROBLEM_H_
#define EXPENSIVESAMPLINGPROBLEM_H_

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ExpensiveSamplingProblem : public SamplingProblem {
    public:

      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target);

      ~ExpensiveSamplingProblem();
      
    private:
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
