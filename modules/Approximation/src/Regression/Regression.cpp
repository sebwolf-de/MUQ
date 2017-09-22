#include "MUQ/Approximation/Regression/Regression.h"

using namespace muq::Modeling;
using namespace muq::Approximation;

Regression::Regression() : WorkPiece() {
  // initalize the algebra
  algebra = std::make_shared<AnyAlgebra>();
}

void Regression::EvaluateImpl(ref_vector<boost::any> const& inputs) {}

//void Regression::Fit(ref_vector<boost::any> const& xs, ref_vector<boost::any> const& ys, boost::any center) const {}

//void Regression::Fit(ref_vector<boost::any> const& xs, ref_vector<boost::any> const& ys) const {}
