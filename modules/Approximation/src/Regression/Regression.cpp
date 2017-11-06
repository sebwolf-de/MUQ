#include "MUQ/Approximation/Regression/Regression.h"

#include "MUQ/Approximation/Regression/Monomial.h"
#include "MUQ/Approximation/Regression/Hermite.h"
#include "MUQ/Approximation/Regression/Legendre.h"

using namespace muq::Modeling;
using namespace muq::Approximation;

Regression::Regression(unsigned int const order, Regression::PolynomialBasis const& basis) : WorkPiece(), order(order) {
  // initalize the algebra
  algebra = std::make_shared<AnyAlgebra>();

  // initalize the polynomial basis
  switch( basis ) {
  case PolynomialBasis::MonomialBasis: {
    poly = std::make_shared<Monomial>();
    break;
  }
  case PolynomialBasis::HermiteBasis: {
    poly = std::make_shared<Hermite>();
    break;
  }
  default: {
    poly = std::make_shared<Legendre>();
    break;
  }
  }
}

void Regression::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // if there are no points ... just return with an empty outputs
  if(inputs.size()==0) { return; }

  // make sure we can compute the Vandermonde matrix
  assert(multi);

  std::vector<boost::any> centered(inputs.size());
  assert(currentRadius>0.0);
  const boost::any normalize = 1.0/currentRadius;
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    // center
    centered[i] = algebra->Subtract(inputs[i].get(), currentCenter);

    // normalize
    centered[i] = algebra->Multiply(normalize, centered[i]);
  }

  // get the Vandermonde matrix of the inputs
  const Eigen::MatrixXd vand = VandermondeMatrix(centered);
  assert(coeff.cols()==vand.cols());

  // compute the regression polynomial
  outputs.resize(1);
  outputs[0] = (Eigen::MatrixXd)(coeff*vand.transpose());
}

int Regression::NumInterpolationPoints() const {
  if( multi ) {
    return multi->Size();
  }

  std::cerr << std::endl << std::endl << "ERROR: Not able to compute the number of points required for interpolation" <<
    std::endl << "\tPolynomialRegressor.cpp NumInterpolationPoints()" << std::endl;
  assert(false);
  
  return -1;
}
