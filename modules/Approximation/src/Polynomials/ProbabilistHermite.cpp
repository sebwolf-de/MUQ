#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"

using namespace muq::Approximation;


double ProbabilistHermite::alpha(unsigned int k, double x) const {
  return -x;
}

double ProbabilistHermite::beta(unsigned int k, double x) const {
  return k;
}

double ProbabilistHermite::phi0(double x) const {
  return 1.0;
}

double ProbabilistHermite::phi1(double x) const {
  return x;
}

REGISTER_POLYNOMIAL_FAMILY(ProbabilistHermite)
