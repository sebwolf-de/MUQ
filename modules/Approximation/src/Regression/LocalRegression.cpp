#include "MUQ/Approximation/Regression/LocalRegression.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

LocalRegression::LocalRegression(std::shared_ptr<WorkPiece> function, pt::ptree const& pt) : WorkPiece(1, std::vector<std::string>(1, typeid(Eigen::VectorXd).name())), kn(pt.get<unsigned int>("LocalRegression.NumNeighbors")) { // can only have one input and output
  // create a cache of model evaluations
  cache = std::make_shared<FlannCache>(function);

  // create a regression object
  reg = std::make_shared<Regression>(pt.get<unsigned int>("LocalRegression.Order", 2));
}

LocalRegression::~LocalRegression() {}

void LocalRegression::FitRegression(boost::any const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  std::vector<Eigen::VectorXd> result;
  cache->NearestNeighbors(input, kn, neighbors, result);

  // fit the regression
  reg->Fit(neighbors, result);
}

void LocalRegression::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // fit the regressor
  FitRegression(inputs[0]);

  // evaluate the regressor
  outputs.resize(1);
  outputs[0] = (Eigen::VectorXd)boost::any_cast<Eigen::MatrixXd const&>(reg->Evaluate(inputs) [0]).col(0);
}

unsigned int LocalRegression::CacheSize() const {
  assert(cache);
  return cache->Size();
}
