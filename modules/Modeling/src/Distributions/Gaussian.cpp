#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(boost::any const& obj, Gaussian::Mode const mode) :
  Distribution(),
  mode(mode),
  dim(algebra->Size(obj,0)),
  mean(mode==Gaussian::Mode::Mean? obj : boost::none),
  cov(mode==Gaussian::Mode::Covariance? SaveCovPrec(obj) : boost::none),
  prec(mode==Gaussian::Mode::Precision? SaveCovPrec(obj) : boost::none) {}

Gaussian::Gaussian(boost::any const& mean, boost::any const& obj, Gaussian::Mode const mode) :
  Distribution(),
  mode(mode),
  dim(algebra->Size(obj,0)),
  mean(mean),
  cov(mode==Gaussian::Mode::Covariance? SaveCovPrec(obj) : boost::none),
  prec(mode==Gaussian::Mode::Precision? SaveCovPrec(obj) : boost::none) {}

Gaussian::~Gaussian() {}

boost::any Gaussian::SaveCovPrec(boost::any const& in) {
  // for Eigen matrices, we can pre compute the cholesky
  if( in.type()==typeid(Eigen::MatrixXd) ) {
    Eigen::LLT<Eigen::MatrixXd> chol;
    chol.compute(boost::any_cast<Eigen::MatrixXd const&>(in));
    
    return chol;
  }

  return in;
}

double Gaussian::LogDensityImpl(ref_vector<boost::any> const& inputs) const {
  boost::any delta = inputs[0].get();
  if( mean ) { delta = algebra->Subtract(delta, *mean); }

  switch( mode ) {
  case Gaussian::Mode::Covariance: {
    return -0.5*algebra->InnerProduct(delta, algebra->ApplyInverse(*cov, delta));
  }
  case Gaussian::Mode::Precision: {
    return -0.5*algebra->InnerProduct(delta, algebra->Apply(*prec, delta));
  }
  case Gaussian::Mode::Mean: {
    return -0.5*algebra->InnerProduct(delta, delta);
  }
  default: {
    assert(false);
    return 0.0;
  }
  }
}

boost::any Gaussian::SampleImpl(ref_vector<boost::any> const& inputs) {
  // make sure the dimension is nonzero
  assert(dim>0);

  const Eigen::VectorXd stdnrm = RandomGenerator::GetNormal(dim);

  switch( mode ) {
  case Gaussian::Mode::Covariance: {
    if( !covSqrt ) { covSqrt = algebra->SquareRoot(*cov); }
    
    // if one dimensional, return a double
    if( dim==1 ) {
      assert(covSqrt->type()==typeid(double));
      return stdnrm(0)*boost::any_cast<double>(*covSqrt);
    }

    // otherwise return an Eigen::VectorXd
    return mean? algebra->Add(algebra->Apply(*covSqrt, stdnrm), *mean) : algebra->Apply(*covSqrt, stdnrm);
  }
  case Gaussian::Mode::Precision: {
    if( !precSqrt ) { precSqrt = algebra->SquareRoot(*prec); }

    // if one dimensional, return a double
    if( dim==1 ) {
      assert(precSqrt->type()==typeid(double));
      return stdnrm(0)/boost::any_cast<double>(*precSqrt);
    }
    
    // otherwise return an Eigen::VectorXd
    return mean? algebra->Add(algebra->ApplyInverse(*precSqrt, stdnrm), *mean) : algebra->ApplyInverse(*precSqrt, stdnrm);
  }
  case Gaussian::Mode::Mean: {
    if( dim==1 ) {
      assert(mean->type()==typeid(double));
      return stdnrm(0)+boost::any_cast<double>(*mean) ;
    }
    return algebra->Add(stdnrm, *mean);
  }
  default: {
    assert(false);
    return boost::none;
  }
  }
}

unsigned int Gaussian::Dimension() const {
  return dim;
}
