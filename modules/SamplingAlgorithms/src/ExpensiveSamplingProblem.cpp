#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, pt::ptree& pt) : SamplingProblem(target) {
  // can only have one input
  assert(target->numInputs==1);
  
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));

  beta = std::pair<double, double>(pt.get<double>("BetaScale", 0.0), -pt.get<double>("BetaExponent", RAND_MAX));
  assert(beta.second<0.0);

  phi = pt.get<double>("StructuralScaling", 1.0);

  lambda = pt.get<double>("PoisednessConstant", 10.0);
  
  gamma = std::pair<double, double>(pt.get<double>("GammaScale", 1.0), pt.get<double>("GammaExponent", 1.0));
  assert(gamma.first>0.0);
  assert(gamma.second>0.0);
}

double ExpensiveSamplingProblem::LogDensity(unsigned int const step, std::shared_ptr<SamplingState> state) {
  std::vector<Eigen::VectorXd> neighbors, results;
  RefineSurrogate(step, state, neighbors, results);
  
  return reg->EvaluateRegressor(state->state[0], neighbors, results) (0);
}

void ExpensiveSamplingProblem::RefineSurrogate(unsigned int const step, std::shared_ptr<SamplingState> state, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) {
  while( reg->CacheSize()<reg->kn ) {
    auto gauss = std::make_shared<muq::Modeling::Gaussian>(state->state[0]);
    reg->Add(gauss->Sample());
  }

  // get the nearest neighbors
  reg->NearestNeighbors(state->state[0], neighbors, results);
  assert(neighbors.size()==results.size());

  // get the error indicator
  std::pair<Eigen::VectorXd, double> error = reg->ErrorIndicator(state->state[0], neighbors);
  
  // BETA refinement
  if( RandomGenerator::GetUniform()<beta.first*std::pow((double)step, beta.second) ) {
    RefineSurrogate(error.first, neighbors, results);
    return;
  }
  
  // check to see if we should increment the level
  if( step>phi*std::pow((double)level, 2.0*gamma.second) ) { ++level; }

  if( error.second>lambda*std::sqrt((double)reg->kn)*gamma.first*std::pow((double)level, -gamma.second) ) {
    RefineSurrogate(error.first, neighbors, results);
    return;
  }
}

void ExpensiveSamplingProblem::RefineSurrogate(Eigen::VectorXd const& point, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) const {
  const Eigen::VectorXd& result = reg->Add(point);
  neighbors[neighbors.size()-1] = point;
  results[results.size()-1] = result;
}

unsigned int ExpensiveSamplingProblem::CacheSize() const { return reg->CacheSize(); }
