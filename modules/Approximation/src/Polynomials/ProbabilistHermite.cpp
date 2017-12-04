#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"

using namespace muq::Approximation;


double ProbabilistHermite::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;
    
    double c = 1.0;
    for(int k=polyOrder; k>polyOrder-derivOrder; --k)
        c *= k;
    
    return c*PolynomialEvaluate(polyOrder-derivOrder, x);
    
}

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
