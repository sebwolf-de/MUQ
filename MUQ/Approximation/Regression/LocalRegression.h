#ifndef LOCALREGRESSION_H_
#define LOCALREGRESSION_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Approximation {
    class LocalRegression : public muq::Modeling::WorkPiece {
    public:

      LocalRegression();

      ~LocalRegression();
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Approximation
} // namespace muq


#endif
