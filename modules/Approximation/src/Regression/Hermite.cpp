#include "MUQ/Approximation/Regression/Hermite.h"

using namespace muq::Approximation;

Hermite::Hermite() : Polynomial() {}

Hermite::~Hermite() {}

double Hermite::alpha(unsigned int k, double x) const {
  return -2.0*x;
}

double Hermite::beta(unsigned int k, double x) const {
  return 2.0*k;
}

double Hermite::phi0(double x) const {
  return 1.0;
}

double Hermite::phi1(double x) const {
  return 2.0*x;
}
