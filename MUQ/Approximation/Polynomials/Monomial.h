#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

namespace muq {
  namespace Approximation {
    /// A monomial polynomial (\f$1\f$, \f$x\f$, \f$x^2\f$, ect. ...)
    /**
       This is a simple polynomial basis but could cause conditioning problems in some cases ...
     */
    class Monomial : public IndexedScalarBasis {
    public:

      Monomial();

      virtual ~Monomial();

      virtual double BasisEvaluate(int const order, double const x) const override;

      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

    };
  } // namespace Approximation
} // namespace muq

#endif
