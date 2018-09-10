#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, pt::ptree pt) : SamplingProblem(target) {
  SetUp(pt);

  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));
}

#if MUQ_HAS_PARCER
ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> comm) : SamplingProblem(target) {
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

  tau0 = pt.get<double>("FirstLevelLength", 1.0);

  gamma = std::pair<double, double>(pt.get<double>("GammaScale", 1.0), pt.get<double>("GammaExponent", 1.0));
  assert(gamma.first>0.0);
  assert(gamma.second>0.0);

  nu = pt.get<double>("TargetMax", 1.0);

  eta = std::pair<double, double>(pt.get<double>("EtaScale", 1.0), pt.get<double>("EtaExponent", 1.0));
}

double ExpensiveSamplingProblem::LogDensity(unsigned int const step, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) {
  std::vector<Eigen::VectorXd> neighbors, results;
  RefineSurrogate(step, state, neighbors, results);
  assert(neighbors.size()==results.size());
  assert(neighbors.size()>=reg->kn);
  /*if( type==AbstractSamplingProblem::SampleType::Accepted ) {
    // refine and get nearest neighbors
    RefineSurrogate(step, state, neighbors, results);
  } else {
    // get nearest neighbors
    reg->NearestNeighbors(state->state[0], neighbors, results);
    }*/

  // set cumulative refinement
  state->meta["cumulative kappa refinement"] = cumkappa;
  state->meta["cumulative beta refinement"] = cumbeta;
  state->meta["cumulative gamma refinement"] = cumgamma;
#if !MUQ_HAS_MPI
  assert(cumbeta+cumgamma+cumkappa==reg->CacheSize());
#endif

  return reg->EvaluateRegressor(state->state[0],
                                std::vector<Eigen::VectorXd>(neighbors.begin(), neighbors.begin()+reg->kn),
                                std::vector<Eigen::VectorXd>(results.begin(), results.begin()+reg->kn)) (0);
}

void ExpensiveSamplingProblem::CheckNumNeighbors(std::shared_ptr<SamplingState> state) {
  while( reg->CacheSize()<reg->kn ) {
    auto gauss = std::make_shared<muq::Modeling::Gaussian>(state->state[0]);
    const Eigen::VectorXd& x = gauss->Sample();
    reg->Add(x);
    UpdateGlobalRadius();
    ++cumkappa;
  }
}

void ExpensiveSamplingProblem::CheckNeighbors(
  std::shared_ptr<SamplingState> state,
  std::vector<Eigen::VectorXd>& neighbors,
  std::vector<Eigen::VectorXd>& results) const {
    // check to see if this state has nearenst neighbors
    const auto& check_neighbors = state->meta.find("nearest neighbors");
    const auto& check_results = state->meta.find("nearest neighbor results");

    if( check_neighbors!=state->meta.end() && check_results!=state->meta.end() ) { // if it does, get them
      neighbors = boost::any_cast<std::vector<Eigen::VectorXd>&>(check_neighbors->second);
      results = boost::any_cast<std::vector<Eigen::VectorXd>&>(check_results->second);
      assert(neighbors.size()==results.size());
    } else { // if it does not, find them
      // get the nearest neighbors
      reg->NearestNeighbors(state->state[0], neighbors, results);
      assert(neighbors.size()==results.size());
      state->meta["nearest neighbors"] = neighbors;
      state->meta["nearest neighbor results"] = results;
    }
}

double ExpensiveSamplingProblem::ErrorThreshold(
  unsigned int const step,
  double const radius,
  double const approxLogTarg) const {
    // a small neggest so rho is always less than one
    const double eps = 1.0e-12;

    // the tail scaling constant
    const double rho = std::fmin(1.0-eps, std::exp(approxLogTarg)/nu);

    // the tail scaling exponent
    const double etaexp = eta.first*std::pow(radius/radius_max, eta.second);

    // compute the error threshold
    return std::pow(rho, etaexp)*gamma.first*std::pow((double)level, -gamma.second);
}

void ExpensiveSamplingProblem::RefineSurrogate(
  std::shared_ptr<SamplingState> state,
  double const radius,
  std::vector<Eigen::VectorXd>& neighbors,
  std::vector<Eigen::VectorXd>& results) {
    // compute the poisedness constant
    const std::tuple<Eigen::VectorXd, double, unsigned int>& lambda
        = reg->PoisednessConstant(state->state[0], neighbors);

    // if the point is already in the cache
    if( reg->InCache(std::get<0>(lambda)) ) {
      // choose a random point
      Eigen::VectorXd point = Eigen::VectorXd::Random(state->state[0].size());
      point *= RandomGenerator::GetUniform()*radius/point.norm();
      point += state->state[0];

      // find the farthest point
      int index = 0;
      double dist = 0.0;
      for( unsigned int i=0; i<reg->kn; ++i ) {
        double newdist = (neighbors[i]-state->state[0]).norm();
        if( newdist>dist ) { dist = newdist; index = i; }
      }

      assert(dist>0.0);
      assert(dist<=radius); // max is the diameter

      RefineSurrogate(point, index, neighbors, results);
    } else {
      RefineSurrogate(std::get<0>(lambda), std::get<2>(lambda), neighbors, results);
    }

  state->meta["nearest neighbors"] = neighbors;
  state->meta["nearest neighbor results"] = results;
}

void ExpensiveSamplingProblem::RefineSurrogate(
  unsigned int const step,
  std::shared_ptr<SamplingState> state,
  std::vector<Eigen::VectorXd>& neighbors,
  std::vector<Eigen::VectorXd>& results) {
  // make sure we have enough points
  CheckNumNeighbors(state);

  // check to see if this state has nearenst neighbors
  CheckNeighbors(state, neighbors, results);

  // get the error indicator (using the first kn neighbors)
  double error, radius;
  std::tie(error, radius) = reg->ErrorIndicator(state->state[0], std::vector<Eigen::VectorXd>(neighbors.begin(), neighbors.begin()+reg->kn));

  // compute the error threshold
  const double approxLogTarg = reg->EvaluateRegressor(state->state[0], neighbors, results) (0);
  const double threshold = ErrorThreshold(step, (state->state[0]-reg->CacheCentroid()).norm(), approxLogTarg);
  state->meta["error indicator"] = error;
  state->meta["error threshold"] = threshold;

  // random refinement
  if( RandomGenerator::GetUniform()<beta.first*std::pow((double)step, beta.second) ) {
    RefineSurrogate(state, radius, neighbors, results);
    ++cumbeta;

    return;
  }

  // check to see if we should increment the level
  if( step>tau0*std::pow((double)level, 2.0*gamma.second) ) { ++level; }

  if( error>threshold ) {
    RefineSurrogate(state, radius, neighbors, results);
    ++cumgamma;
    return;
  }
}

void ExpensiveSamplingProblem::RefineSurrogate(Eigen::VectorXd const& point, unsigned int const index, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) {
  assert(!reg->InCache(point));
  const Eigen::VectorXd& result = reg->Add(point);

  assert(neighbors.size()==results.size());
  assert(index<neighbors.size());

  neighbors.push_back(point);
  results.push_back(result);
  std::iter_swap(neighbors.end()-1, neighbors.begin()+index);
  std::iter_swap(results.end()-1, results.begin()+index);

  UpdateGlobalRadius();
}

void ExpensiveSamplingProblem::UpdateGlobalRadius() {
  radius_max = 0.0;
  for( unsigned int i=0; i<CacheSize(); ++i ) {
    const double r = (reg->CacheCentroid()-reg->CachePoint(i)).norm();
    radius_max = std::max(radius_max, r);
  }
}

unsigned int ExpensiveSamplingProblem::CacheSize() const { return reg->CacheSize(); }
