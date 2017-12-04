#ifndef LAGUERRE_H_
#define LAGUERRE_H_

#include "MUQ/Approximation/Polynomials/Polynomial.h"

namespace muq {
    namespace Approximation {
        class Laguerre : public Polynomial {
        public:
            
            /**
               Laguerre polynomials are orthogonal over \f$[0,\infty)\f$ with respect to \f$\exp(-x)x^a\f$.
            */
            Laguerre(const double aIn=0.0) : a(aIn){};
            
            virtual ~Laguerre() = default;
            
            virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;
            
            virtual double Normalization(unsigned int polyOrder) const override;
            
        private:
            
            const double a;
            
            virtual double alpha(unsigned int k, double x) const override;
            virtual double beta(unsigned int k, double x) const override;
            
            virtual double phi0(double x) const override;
            virtual double phi1(double x) const override;
        };
    } // namespace Approximation
} // namespace muq

#endif
