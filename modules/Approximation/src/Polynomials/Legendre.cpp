#include "MUQ/Approximation/Polynomials/Legendre.h"

using namespace muq::Approximation;

Legendre::Legendre() : Polynomial() {}

Legendre::~Legendre() {}

double Legendre::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    if(derivOrder==1){
        return polyOrder / (x * x - 1.0) * (x * PolynomialEvaluate(polyOrder, x) - PolynomialEvaluate(polyOrder - 1, x));
    }else{
        // Use the fact that dp_{n+1}/dx = (2n+1) p_n + dp_{n-1} /dx
        return (2*(polyOrder-1) + 1) * DerivativeEvaluate(polyOrder-1, derivOrder-1, x) + DerivativeEvaluate(polyOrder-2, derivOrder, x);
    }
}


double Legendre::alpha(unsigned int k, double x) const {
  return -(2.0*(double)k+1.0)*x/((double)k+1.0);
}

double Legendre::beta(unsigned int k, double x) const {
  return (double)k/((double)k+1.0);
}

double Legendre::phi0(double x) const {
  return 1.0;
}

double Legendre::phi1(double x) const {
  return x;
}

double Legendre::Normalization(unsigned int polyOrder) const {
    return 2.0/(2.0*polyOrder + 1.0);
}


REGISTER_POLYNOMIAL_FAMILY(Legendre)
