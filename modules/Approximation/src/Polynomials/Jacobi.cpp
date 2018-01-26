#include "MUQ/Approximation/Polynomials/Jacobi.h"

using namespace muq::Approximation;


double Jacobi::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {
   
    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = std::tgamma(a+b+polyOrder+1+derivOrder) / (std::pow(2.0, derivOrder) * std::tgamma(a+b+polyOrder+1));
    
    return c * Jacobi(a + derivOrder, b+derivOrder).PolynomialEvaluate(polyOrder-derivOrder, x);
}


double Jacobi::alpha(unsigned int polyOrder, double x) const {

    if(polyOrder==0)
        return -1.0*(0.5*(a+b)+1)*x - 0.5*(a-b);
    
    const double den1 = 2.0*(polyOrder+1)*(polyOrder + a + b  + 1.0);
    const double den2 = (2.0*polyOrder + a + b) * den1;
    
    const double An = (2.0*polyOrder + a + b + 1.0)*(2.0*polyOrder + a + b + 2.0) / den1;
    const double Bn = (a*a-b*b) * (2.0*polyOrder + a + b + 1.0) / den2;
        
    return -An*x - Bn;
}

double Jacobi::beta(unsigned int polyOrder, double x) const {

    const double den = (polyOrder+1)*(polyOrder + a + b  + 1.0)*(2.0*polyOrder + a + b);

    const double Cn = (polyOrder + a)*(polyOrder + b)*(2.0*polyOrder + a + b + 2.0) / den;
        
    return Cn;
}

double Jacobi::phi0(double x) const {
    return 1.0;
}

double Jacobi::phi1(double x) const {
    return (0.5*(a+b)+1.0)*x + 0.5*(a-b);
}

double Jacobi::Normalization(unsigned int polyOrder) const {
    return (std::pow(2.0, a+b+1) / (2*polyOrder + a + b +1)) * std::tgamma(polyOrder+a+1)*std::tgamma(polyOrder+b+1)/(std::tgamma(polyOrder+a+b+1) * std::tgamma(polyOrder + 1));
}


REGISTER_POLYNOMIAL_FAMILY(Jacobi)
