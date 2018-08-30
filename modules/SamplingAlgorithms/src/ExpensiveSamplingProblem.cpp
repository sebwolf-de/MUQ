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
  SetUp(pt);

  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));
}

#if MUQ_HAS_PARCER
ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree& pt, std::shared_ptr<parcer::Communicator> comm) : SamplingProblem(target) {
  SetUp(pt);

  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")), comm);
}
#endif

void ExpensiveSamplingProblem::SetUp(boost::property_tree::ptree& pt) {
  // can only have one input
  assert(target->numInputs==1);

  beta = std::pair<double, double>(pt.get<double>("BetaScale", 0.0), -pt.get<double>("BetaExponent", RAND_MAX));
  assert(beta.second<0.0);

  phi = pt.get<double>("StructuralScaling", 1.0);

  lambda = pt.get<double>("PoisednessConstant", 10.0);

  gamma = std::pair<double, double>(pt.get<double>("GammaScale", 1.0), pt.get<double>("GammaExponent", 1.0));
  assert(gamma.first>0.0);
  assert(gamma.second>0.0);
}

double ExpensiveSamplingProblem::LogDensity(unsigned int const step, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) {
  std::vector<Eigen::VectorXd> neighbors, results;
  RefineSurrogate(step, state, neighbors, results);
  /*if( type==AbstractSamplingProblem::SampleType::Accepted ) {
    // refine and get nearest neighbors
    RefineSurrogate(step, state, neighbors, results);
  } else {
    // get nearest neighbors
    reg->NearestNeighbors(state->state[0], neighbors, results);
    }*/

  // set cumulative refinement
  state->meta["cumulative beta refinement"] = cumbeta;
  state->meta["cumulative gamma refinement"] = cumgamma;
  state->meta["cumulative kappa refinement"] = cumkappa;
#if !MUQ_HAS_MPI
  assert(cumbeta+cumgamma+cumkappa==reg->CacheSize());
#endif

  return reg->EvaluateRegressor(state->state[0], neighbors, results) (0);
}

void ExpensiveSamplingProblem::RefineSurrogate(unsigned int const step, std::shared_ptr<SamplingState> state, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) {
  while( reg->CacheSize()<reg->kn ) {
    auto gauss = std::make_shared<muq::Modeling::Gaussian>(state->state[0]);
    reg->Add(gauss->Sample());
    ++cumkappa;
  }

  // get the nearest neighbors
  reg->NearestNeighbors(state->state[0], neighbors, results);
  assert(neighbors.size()==results.size());

  // get the error indicator
  const std::tuple<Eigen::VectorXd, double, unsigned int>& error = reg->ErrorIndicator(state->state[0], neighbors);
  state->meta["error indicator"] = std::get<1>(error);
  state->meta["error threshold"] = std::numeric_limits<double>::quiet_NaN();

  // BETA refinement
  if( RandomGenerator::GetUniform()<beta.first*std::pow((double)step, beta.second) ) {
    const std::tuple<Eigen::VectorXd, double, unsigned int>& lambda = reg->PoisednessConstant(state->state[0], neighbors);
    RefineSurrogate(std::get<0>(lambda), std::get<2>(lambda), neighbors, results);
    ++cumbeta;
    return;
  }

  // check to see if we should increment the level
  if( step>phi*std::pow((double)level, 2.0*gamma.second) ) { ++level; }

  const double threshold = lambda*std::sqrt((double)reg->kn)*gamma.first*std::pow((double)level, -gamma.second);
  state->meta["error threshold"] = threshold;
  if( std::get<1>(error)>threshold ) {
    const std::tuple<Eigen::VectorXd, double, unsigned int>& lambda = reg->PoisednessConstant(state->state[0], neighbors);
    RefineSurrogate(std::get<0>(lambda), std::get<2>(lambda), neighbors, results);
    ++cumgamma;
    return;
  }
}

void ExpensiveSamplingProblem::RefineSurrogate(Eigen::VectorXd const& point, unsigned int const index, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) const {
  const Eigen::VectorXd& result = reg->Add(point);
  neighbors[index] = point;
  results[index] = result;
}

unsigned int ExpensiveSamplingProblem::CacheSize() const { return reg->CacheSize(); }
