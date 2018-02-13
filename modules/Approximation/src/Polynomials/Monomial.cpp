#include "MUQ/Approximation/Polynomials/Monomial.h"

using namespace muq::Approximation;

Monomial::Monomial() : Polynomial() {}

Monomial::~Monomial() {}

double Monomial::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = 1.0;
    for(int k=polyOrder; k>polyOrder-derivOrder; --k)
        c *= k;

    return c*std::pow(x, polyOrder-derivOrder);

}

double Monomial::BasisEvaluate(int const order, double const x) const {
    return std::pow(x, order);
}

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

REGISTER_SCALARBASIS_FAMILY(Monomial)
