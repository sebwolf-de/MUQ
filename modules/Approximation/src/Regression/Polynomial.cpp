#include "MUQ/Approximation/Regression/Polynomial.h"

using namespace muq::Modeling;
using namespace muq::Approximation;

Polynomial::Polynomial() :
  WorkPiece(std::vector<std::string>({typeid(unsigned int).name(), typeid(double).name()}), // input types (order, point)
	    std::vector<std::string>({typeid(double).name()})) // output times (polynomial evaluation)
{}

Polynomial::~Polynomial() {}

void Polynomial::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // extract the inputs
  const unsigned int order = boost::any_cast<unsigned int>(inputs[0]);
  const double x = boost::any_cast<double>(inputs[1]);

  // evaluate the polynomial
  outputs.resize(1);
  outputs[0] = PolynomialEvaluate(order, x);
}

double Polynomial::PolynomialEvaluate(int const order, double const x) const {
  double bkp2 = 0.0;
  double bkp1 = 0.0;
  double bk = 1.0;

  for( int k=order-1; k>=0; k-- ) {
    // increment
    bkp2 = bkp1;
    bkp1 = bk;

    // compute new bk
    bk = -alpha(k, x)*bkp1 - beta(k+1, x)*bkp2;
  }

  return bk*phi0(x) + bkp1*(phi1(x)+alpha(0, x)*phi0(x));
}
