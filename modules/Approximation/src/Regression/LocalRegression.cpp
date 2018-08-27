#include "MUQ/Approximation/Regression/LocalRegression.h"

#include <parcer/Eigen.h>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

LocalRegression::LocalRegression(std::shared_ptr<ModPiece> function, pt::ptree& pt) : ModPiece(function->inputSizes, function->outputSizes), kn(pt.get<unsigned int>("NumNeighbors")) { 
  SetUp(function, pt);
}

#if MUQ_HAS_PARCER
LocalRegression::LocalRegression(std::shared_ptr<muq::Modeling::ModPiece> function, boost::property_tree::ptree& pt, std::shared_ptr<parcer::Communicator> comm) : ModPiece(function->inputSizes, function->outputSizes), kn(pt.get<unsigned int>("NumNeighbors")), comm(comm) {
  SetUp(function, pt);
}
#endif

void LocalRegression::SetUp(std::shared_ptr<muq::Modeling::ModPiece> function, boost::property_tree::ptree& pt) {
  // can only have one input and output
  assert(inputSizes.size()==1);
  assert(outputSizes.size()==1);
  
  // create a cache of model evaluations
  cache = std::make_shared<FlannCache>(function);
  
  // create a regression object
  pt.put<std::string>("PolynomialBasis", pt.get<std::string>("PolynomialBasis", "Legendre")); // set default to Legendre
  pt.put<unsigned int>("Order", pt.get<unsigned int>("Order", 2)); // set default order to 2
  pt.put<unsigned int>("InputSize", function->inputSizes(0));
  reg = std::make_shared<Regression>(pt);
}

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

Eigen::VectorXd LocalRegression::Add(Eigen::VectorXd const& input) const {
  assert(cache);
  const Eigen::VectorXd& result = cache->Add(input);
  
#if MUQ_HAS_PARCER
  if( comm ) {
    for( unsigned int i=0; i<comm->GetSize(); ++i ) {
      if( i==comm->GetRank() ) { continue; }
      
      parcer::SendRequest sendReq;
      comm->Isend(std::pair<Eigen::VectorXd, Eigen::VectorXd>(input, result), i, tagSingle, sendReq);
    }

    Probe();
  }
#endif  
  
  assert(cache);
  return result;
}

void LocalRegression::Add(std::vector<Eigen::VectorXd> const& inputs) const {
  assert(cache);
  const std::vector<Eigen::VectorXd>& results = cache->Add(inputs);

#if MUQ_HAS_PARCER
  if( comm ) {
    for( unsigned int i=0; i<comm->GetSize(); ++i ) {
      if( i==comm->GetRank() ) { continue; }
      
      parcer::SendRequest sendReq;
      comm->Isend(std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd> >(inputs, results), i, tagMulti, sendReq);
    }
  }
#endif
}

std::tuple<Eigen::VectorXd, double, unsigned int> LocalRegression::PoisednessConstant(Eigen::VectorXd const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  cache->NearestNeighbors(input, kn, neighbors);

  return PoisednessConstant(input, neighbors);
}

std::tuple<Eigen::VectorXd, double, unsigned int> LocalRegression::PoisednessConstant(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors) const {
  assert(reg);
  std::pair<Eigen::VectorXd, double> lambda = reg->PoisednessConstant(neighbors, input);

  double dist = RAND_MAX;
  unsigned int index = 0;
  for( unsigned int i=0; i<neighbors.size(); ++i ) {
    const double d = (lambda.first-neighbors[i]).norm();
    if( d<dist ) { dist = d; index = i; }
  }

  return std::tuple<Eigen::VectorXd, double, unsigned int>(lambda.first, lambda.second, index);
}

std::tuple<Eigen::VectorXd, double, unsigned int> LocalRegression::ErrorIndicator(Eigen::VectorXd const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  cache->NearestNeighbors(input, kn, neighbors);

  return ErrorIndicator(input, neighbors);
}

std::tuple<Eigen::VectorXd, double, unsigned int> LocalRegression::ErrorIndicator(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors) const {
  // get the poisedness constant
  std::tuple<Eigen::VectorXd, double, unsigned int> lambda = PoisednessConstant(input, neighbors);

  // update the error indicator
  std::get<1>(lambda) *= std::sqrt((double)kn)*std::pow((*(neighbors.end()-1)-input).norm(), (double)reg->order+1.0);

  return lambda;
}

void LocalRegression::NearestNeighbors(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd>& neighbors) const {
  assert(cache);
  cache->NearestNeighbors(input, kn, neighbors);
}

void LocalRegression::NearestNeighbors(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& result) const {
  assert(cache);
  cache->NearestNeighbors(input, kn, neighbors, result);
}

Eigen::VectorXd LocalRegression::EvaluateRegressor(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors, std::vector<Eigen::VectorXd> const& result) const {
  // fit the regression
  reg->Fit(neighbors, result, input);

  // evaluate the regressor
  return (Eigen::VectorXd)boost::any_cast<Eigen::MatrixXd const&>(reg->Evaluate(input) [0]).col(0);  
}

#if MUQ_HAS_PARCER
void LocalRegression::Probe() const {
  if( !comm ) { return; }

  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    if( i==comm->GetRank() ) { continue; }

    { // get single adds
      parcer::RecvRequest recvReq;
      bool has = comm->Iprobe(i, tagSingle, recvReq);
      while( has ) {
	const std::pair<Eigen::VectorXd, Eigen::VectorXd>& other = comm->Recv<std::pair<Eigen::VectorXd, Eigen::VectorXd> >(i, tagSingle);
	cache->Add(other.first, other.second);
	has = comm->Iprobe(i, tagSingle, recvReq);
      }
    }

    { // get multi adds
      parcer::RecvRequest recvReq;
      bool has = comm->Iprobe(i, tagMulti, recvReq);
      while( has ) {
	const std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>>& other = comm->Recv<std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd> > >(i, tagMulti);
	cache->Add(other.first, other.second);
	has = comm->Iprobe(i, tagMulti, recvReq);
      }
    }
  }
}
#endif
