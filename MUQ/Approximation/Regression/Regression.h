#ifndef REGRESSION_H_
#define REGRESSION_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Approximation {
    class Regression : public muq::Modeling::WorkPiece {
    public:

      Regression();
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
