#include "MUQ/Approximation/Regression/Monomial.h"

using namespace muq::Approximation;

Monomial::Monomial() : Polynomial() {}

Monomial::~Monomial() {}

double Monomial::alpha(unsigned int k, double x) const {
  return -x;
}

double Monomial::beta(unsigned int k, double x) const {
  return 0.0;
}

double Monomial::phi0(double x) const {
  return 1.0;
}

double Monomial::phi1(double x) const {
  return x;
}
