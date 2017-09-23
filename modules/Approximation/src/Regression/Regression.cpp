#include "MUQ/Approximation/Regression/Regression.h"

#include "MUQ/Approximation/Regression/Monomial.h"
#include "MUQ/Approximation/Regression/Hermite.h"
#include "MUQ/Approximation/Regression/Legendre.h"

using namespace muq::Modeling;
using namespace muq::Approximation;

Regression::Regression(unsigned int const dim, unsigned int const order, Regression::PolynomialBasis const& basis) : WorkPiece() {
  // initalize the algebra
  algebra = std::make_shared<AnyAlgebra>();

  // initalize the multi-index
  multi = std::make_shared<MultiIndex>(dim, order);

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

  // get the Vandermonde matrix of the inputs
  const Eigen::MatrixXd vand = VandermondeMatrix(inputs);
  assert(coeff.cols()==vand.cols());

  // compute the regression polynomial
  outputs.resize(1);
  outputs[0] = (Eigen::MatrixXd)(coeff*vand.transpose());
}
