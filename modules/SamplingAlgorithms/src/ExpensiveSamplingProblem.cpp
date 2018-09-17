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
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));

  // must becaused after reg is created
  SetUp(pt);
}

#if MUQ_HAS_PARCER
ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> comm) : SamplingProblem(target) {
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")), comm);

  // must becaused after reg is created
  SetUp(pt);
}
#endif

void ExpensiveSamplingProblem::SetUp(boost::property_tree::ptree& pt) {
  assert(reg);

  // can only have one input
  assert(target->numInputs==1);

  beta = std::pair<double, double>(pt.get<double>("BetaScale", 0.0), -pt.get<double>("BetaExponent", RAND_MAX));
  assert(beta.second<0.0);

  tau0 = pt.get<double>("FirstLevelLength", 1.0);

  gamma = std::pair<double, double>(pt.get<double>("GammaScale", 1.0), pt.get<double>("GammaExponent", 1.0));
  assert(gamma.first>0.0);
  assert(gamma.second>0.0);
  maxGammaRefine = pt.get<unsigned int>("MaximumGammaRefine", reg->kn);

  nu = pt.get<double>("TargetMax", 1.0);

  lambdaScale = pt.get<double>("LambdaScale", 25.0);

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
  state->meta["global radius"] = radius_max;

  return reg->EvaluateRegressor(state->state[0],
                                std::vector<Eigen::VectorXd>(neighbors.begin(), neighbors.begin()+reg->kn),
                                std::vector<Eigen::VectorXd>(results.begin(), results.begin()+reg->kn)) (0);
}

void ExpensiveSamplingProblem::CheckNumNeighbors(std::shared_ptr<SamplingState> state) {
  while( reg->CacheSize()<reg->kn ) {
    auto gauss = std::make_shared<muq::Modeling::Gaussian>(state->state[0], 1.0e-3*Eigen::VectorXd::Ones(state->state[0].size()));
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
    reg->NearestNeighbors(state->state[0], neighbors, results);
    assert(neighbors.size()==results.size());
}

std::tuple<double, double, double> ExpensiveSamplingProblem::ErrorThreshold(
  unsigned int const step,
  double const radius,
  double const approxLogTarg) const {
    // a small neggest so rho is always less than one
    const double eps = 1.0e-12;

    // the tail scaling constant
    const double rho = std::fmin(1.0-eps, std::exp(approxLogTarg-nu));

    // the tail scaling exponent
    const double etaexp = eta.first*std::pow(radius/radius_max, eta.second);

    return std::tuple<double, double, double>(rho, etaexp, gamma.first*std::pow((double)level, -gamma.second));
}

double ExpensiveSamplingProblem::RefineSurrogate(
  std::shared_ptr<SamplingState> state,
  double const radius,
  std::vector<Eigen::VectorXd>& neighbors,
  std::vector<Eigen::VectorXd>& results) {
    // compute the poisedness constant
    const std::tuple<Eigen::VectorXd, double, unsigned int>& lambda = reg->PoisednessConstant(state->state[0], neighbors);

    // if the point is already in the cache
    if( reg->InCache(std::get<0>(lambda)) ) {
      // choose a random point
      Eigen::VectorXd point = Eigen::VectorXd::Random(state->state[0].size());
      point *= RandomGenerator::GetUniform()*radius/point.norm();
      point += state->state[0];

      // find the closest point
      int index = 0;
      double dist = RAND_MAX;
      for( unsigned int i=0; i<reg->kn; ++i ) {
        double newdist = (point-state->state[0]).norm();
        if( newdist<dist ) { dist = newdist; index = i; }
      }

      assert(dist>0.0);
      assert(dist<=radius); // max is the diameter

      RefineSurrogate(point, index, neighbors, results);
    } else {
      RefineSurrogate(std::get<0>(lambda), std::get<2>(lambda), neighbors, results);
    }

  return std::get<1>(lambda);
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

  // random refinement
  if( RandomGenerator::GetUniform()<beta.first*std::pow((double)step, beta.second) ) {
    RefineSurrogate(state, radius, neighbors, results);
    ++cumbeta;

    // recompute the error indicator
    std::tie(error, radius) = reg->ErrorIndicator(state->state[0], std::vector<Eigen::VectorXd>(neighbors.begin(), neighbors.begin()+reg->kn));
  }

  // compute the error threshold
  const double approxLogTarg = reg->EvaluateRegressor(state->state[0], neighbors, results) (0);

  // check to see if we should increment the level
  if( step>tau0*std::pow((double)level, 2.0*gamma.second) ) { ++level; }

  // compute (and store) the error threshold
  double rho, etaexp, threshold;
  std::tie(rho, etaexp, threshold) = ErrorThreshold(step, (state->state[0]-reg->CacheCentroid()).norm(), approxLogTarg);
  threshold *= std::pow(rho, etaexp);
  state->meta["error indicator"] = error;
  state->meta["nearest neighbor ball radius"] = radius;
  state->meta["rho"] = rho;
  state->meta["eta"] = etaexp;
  state->meta["error threshold"] = threshold;

  // structural refinement
  unsigned int n = 0;
  if( error>threshold ) {
    double lambda = RefineSurrogate(state, radius, neighbors, results);
    ++cumgamma; ++n;
    while( lambda>lambdaScale*reg->kn && n<=maxGammaRefine ) {
      lambda = RefineSurrogate(state, radius, neighbors, results);
      ++cumgamma; ++n;
    }
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
    radius_max = std::fmax(radius_max, r);
  }
}

unsigned int ExpensiveSamplingProblem::CacheSize() const { return reg->CacheSize(); }
