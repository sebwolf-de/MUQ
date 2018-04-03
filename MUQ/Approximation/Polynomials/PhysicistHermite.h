#ifndef PHYSICISTHERMITE_H_
#define PHYSICISTHERMITE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

namespace muq {
  namespace Approximation {
    class PhysicistHermite : public OrthogonalPolynomial{
    public:

      /// A Hermite polynomial (\f$1\f$, \f$2x\f$, \f$4x^2-2.0\f$, ect. ...)
      /**
	 Hermite polynomials are orthogonal, which helps with some conditioning problems.   Here we implement the physicists' Hermite polynomials.

	 note: Hermite may be unstable for high orders.
       */
      PhysicistHermite();

      virtual ~PhysicistHermite();

      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

      virtual double Normalization(unsigned int polyOrder) const override;

    private:

      virtual double ak(unsigned int k) const override;
      virtual double bk(unsigned int k) const override;
      virtual double ck(unsigned int k) const override;
      virtual double phi0(double x) const override;
      virtual double phi1(double x) const override;

    };
  } // namespace Approximation
} // namespace muq

#endif
