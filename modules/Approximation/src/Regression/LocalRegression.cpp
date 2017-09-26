#include "MUQ/Approximation/Regression/LocalRegression.h"

using namespace muq::Modeling;
using namespace muq::Approximation;

LocalRegression::LocalRegression(std::shared_ptr<WorkPiece> function) {
  // create a cache of model evaluations
  cache = std::make_shared<FlannCache>(function);
}

LocalRegression::~LocalRegression() {}

void LocalRegression::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  
}
