#ifndef PROBABILISTHERMITE_H_
#define PROBABILISTHERMITE_H_

#include "MUQ/Approximation/Polynomials/Polynomial.h"

namespace muq {
  namespace Approximation {
    class ProbabilistHermite : public Polynomial{
    public:

      /// A Hermite polynomial (\f$1\f$, \f$2x\f$, \f$4x^2-2.0\f$, ect. ...)
      /**
	 Hermite polynomials are orthogonal, which helps with some conditioning problems.   Here we implement the probabilists' Hermite polynomials.

	 note: Hermite may be unstable for high orders.
       */
      ProbabilistHermite() = default;

      virtual ~ProbabilistHermite() = default;

      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;
      
    private:

      /// Implement \f$\alpha_k(x)\f$
      /**
	 @param[in] k The order of the polynomial
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double alpha(unsigned int k, double x) const override;

      /// Implement \f$\beta_k(x)\f$
      /**
	 @param[in] k The order of the polynomial
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double beta(unsigned int k, double x) const override;

      /// Implement \f$\phi_0(x)\f$
      /**
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi0(double x) const override;

      /// Implement \f$\phi_1(x)\f$
      /**
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi1(double x) const override;
    };
  } // namespace Approximation
} // namespace muq

#endif
