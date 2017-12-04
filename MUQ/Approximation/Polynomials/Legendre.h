#ifndef LEGENDRE_H_
#define LEGENDRE_H_

#include "MUQ/Approximation/Polynomials/Polynomial.h"

namespace muq {
  namespace Approximation {
    class Legendre : public Polynomial {
    public:

      /// A Legendre polynomial (\f$1\f$, \f$x\f$, \f$\frac{1}{2}(3x^2-1)\f$, ect. ...)
      /**
	 Legendre polynomials are orthogonal, which helps with some conditioning problems.
       */
      Legendre();

      virtual ~Legendre();
      
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
