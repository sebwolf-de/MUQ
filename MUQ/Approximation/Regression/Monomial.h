#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include "MUQ/Approximation/Regression/Polynomial.h"

namespace muq {
  namespace Approximation {
    /// A monomial polynomial (\f$1\f$, \f$x\f$, \f$x^2\f$, ect. ...)
    class Monomial : public Polynomial {
    public:

      Monomial();

      virtual ~Monomial();
      
    private:
    };
  } // namespace Approximation
} // namespace muq

#endif
