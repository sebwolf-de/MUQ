#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <functional>
#include <string>

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

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
    class Polynomial : public IndexedScalarBasis {
    public:

      /// Create a polynomial
      Polynomial() : IndexedScalarBasis(){};

      virtual ~Polynomial() = default;

      /// Evaluate the specific polynomial type (must be implemented by the child)
      /**
	 Inputs:
	 <ol>
	 <li> The order of the polynomial (unsigned int)
	 <li> The point where we are evaluating the polynomial
	 </ol>
	 \return The polynomial value
       */
      virtual double BasisEvaluate(int const order, double const x) const override;

      /** Evaluate the \f$n^{th}\f$ derivative of a \f$p^{th}\f$ order polynomial.
          @param[in] polyOrder The order \f$p\f$ of the polynomial.
          @param[in] derivOrder The order \f$n\f$ of the derivative.
          @param[in] x The location to evaluate the derivative.
      */
      //virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const = 0;



      /**
         @brief Returns the normalization constant for the polynomial of order \f$p\f$.
         @details Many polynomial families are orthogonal with respect to some measure.  This means that for two polynomials \f$P_n(x)\f$ and \f$P_m(x)\f$ in the same family,
\f[
\int P_n(x) P_m(x) w(x) dx = a_n \delta_{nm},
\f]
where \f$\delta_{mn}\f$ is the Kronecker delta function, \f$w(x)\f$ is a weighting function tied to the polynomial family (e.g., the Askey scheme), and \f$a_n\in\mathbb{R}\f$ is a scalar normalization constant.  This function computes the normalization constant \f$a_n\f$.  Polynomial families that are not orthogonal do not have a normalizing constant \f$a_n\f$ and will throw an exception if this function is called.
         @param[in] polyOrder The order of the polynomial
         @return The normalization constant for a particular polynomial family and its weighting function \f$w(x)\f$.
       */
      virtual double Normalization(unsigned int polyOrder) const;


    private:

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
