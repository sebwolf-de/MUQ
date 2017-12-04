#include "MUQ/Approximation/Polynomials/Legendre.h"

using namespace muq::Approximation;

Legendre::Legendre() : Polynomial() {}

Legendre::~Legendre() {}

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


REGISTER_POLYNOMIAL_FAMILY(Legendre)
