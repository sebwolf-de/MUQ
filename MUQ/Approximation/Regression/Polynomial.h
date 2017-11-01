#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Approximation {
    /// A 1D polynomial (monomial, Hermite, or Legendre)
    /**
       In general, we use an recursive formula to evaluate a \f$d^{th}\f$ degree polynomial using a three term recurrence:
       \f{eqnarray*}{
       p_0(x) &=& \phi_0(x) \\
       p_1(x) &=& \phi_1(x) \\
       p_{k+1}(x) &=& - \alpha_k(x)p_k(x) - \beta_k(x) p_{k-1}(x)
       \f}
       Subclasses specialize for particular polynomials (e.g., Hermite) by implementing \f$\alpha_k(x)\f$, \f$\beta_k(x)\f$, \f$\phi_0(x)\f$, and \f$\phi_1(x)\f$.

       Uses the Clenshaw algorithm from: http://en.wikipedia.org/wiki/Clenshaw_algorithm.
     */
    class Polynomial : public muq::Modeling::WorkPiece {
    public:

      /// Create a polynomial
      Polynomial();

      virtual ~Polynomial();
      
    private:

      /// Evaluate the polynomial at a given point and order
      /**
	 Inputs:
	 <ol>
	 <li> The order of the polynomial (unsigned int)
	 <li> The point where we are evaluating the polynomial (double)
	 </ol>
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs);

      /// Evaluate the specific polynomial type (must be implemented by the child)
      /**
	 Inputs:
	 <ol>
	 <li> The order of the polynomial (unsigned int)
	 <li> The point where we are evaluating the polynomial
	 </ol>
	 \return The polynomial value
       */
      virtual double PolynomialEvaluate(int const order, double const x) const;

      /// Implement \f$\alpha_k(x)\f$
      /**
	 @param[in] k The order of the polynomial
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double alpha(unsigned int k, double x) const = 0;

      /// Implement \f$\beta_k(x)\f$
      /**
	 @param[in] k The order of the polynomial
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double beta(unsigned int k, double x) const = 0;

      /// Implement \f$\phi_0(x)\f$
      /**
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi0(double x) const = 0;

      /// Implement \f$\phi_1(x)\f$
      /**
	 @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi1(double x) const = 0;
    };
  } // namespace Approximation
} // namespace muq

#endif
