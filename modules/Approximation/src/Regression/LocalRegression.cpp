#include "MUQ/Approximation/Regression/LocalRegression.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

LocalRegression::LocalRegression(std::shared_ptr<WorkPiece> function, pt::ptree const& pt) : WorkPiece(1, 1), kn(pt.get<unsigned int>("LocalRegression.NumNeighbors")) { // can only have one input and output
  // create a cache of model evaluations
  cache = std::make_shared<FlannCache>(function);

  // create a regression object
  reg = std::make_shared<Regression>(pt.get<unsigned int>("LocalRegression.Order", 2), "Legendre");
}

LocalRegression::~LocalRegression() {}

void LocalRegression::FitRegression(boost::any const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  std::vector<Eigen::VectorXd> result;
  cache->NearestNeighbors(input, kn, neighbors, result);

  // fit the regression
  reg->Fit(neighbors, result, input);
}

void LocalRegression::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // fit the regressor
  FitRegression(inputs[0]);

  // evaluate the regressor
  const Eigen::VectorXd result = boost::any_cast<Eigen::MatrixXd const&>(reg->Evaluate(inputs) [0]).col(0);

  // store the result as the appropriate output time
  outputs.resize(1);
  switch( result.size() ) {
  case 1: { outputs[0] = result(0); break; }
  case 2: { outputs[0] = (Eigen::Vector2d)result; break; }
  case 3: { outputs[0] = (Eigen::Vector3d)result; break; }
  case 4: { outputs[0] = (Eigen::Vector4d)result; break; }
  default: { outputs[0] = result; }
  }
}

unsigned int LocalRegression::CacheSize() const {
  assert(cache);
  return cache->Size();
}
