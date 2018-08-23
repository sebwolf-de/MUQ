#include "MUQ/Approximation/Regression/LocalRegression.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

LocalRegression::LocalRegression(std::shared_ptr<ModPiece> function, pt::ptree const& pt) : ModPiece(function->inputSizes, function->outputSizes), kn(pt.get<unsigned int>("NumNeighbors")) { // can only have one input and output
  assert(inputSizes.size()==1);
  assert(outputSizes.size()==1);
  
  // create a cache of model evaluations
  cache = std::make_shared<FlannCache>(function);

  // create a regression object
  reg = std::make_shared<Regression>(pt.get<unsigned int>("Order", 2), "Legendre");
}

LocalRegression::~LocalRegression() {}

void LocalRegression::FitRegression(Eigen::VectorXd const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  std::vector<Eigen::VectorXd> result;
  cache->NearestNeighbors(input, kn, neighbors, result);

  // fit the regression
  reg->Fit(neighbors, result, input);
}

void LocalRegression::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // fit the regressor
  FitRegression(inputs[0]);

  // evaluate the regressor
  outputs.resize(1);
  outputs[0] = (Eigen::VectorXd)boost::any_cast<Eigen::MatrixXd const&>(reg->Evaluate(inputs[0].get()) [0]).col(0);
}

unsigned int LocalRegression::CacheSize() const {
  assert(cache);
  return cache->Size();
}

void LocalRegression::Add(std::vector<Eigen::VectorXd> const& inputs) const {
  assert(cache);
  for( auto it : inputs ) { cache->Add(it); }
}
