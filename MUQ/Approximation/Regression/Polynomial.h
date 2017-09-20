#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Approximation {
    /// A 1D polynomial (monomial, Hermite, or Legendre)
    class Polynomial : public muq::Modeling::WorkPiece {
    public:

      Polynomial();

      virtual ~Polynomial();
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs);
      
    };
  } // namespace Approximation
} // namespace muq

#endif
