#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include "MUQ/Approximation/Polynomials/Polynomial.h"

namespace muq {
  namespace Approximation {
    /// A monomial polynomial (\f$1\f$, \f$x\f$, \f$x^2\f$, ect. ...)
    /**
       This is a simple polynomial basis but could cause conditioning problems in some cases ...
     */
    class Monomial : public Polynomial {
    public:

      Monomial();

      virtual ~Monomial();

      virtual double PolynomialEvaluate(int const order, double const x) const override;
      
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
