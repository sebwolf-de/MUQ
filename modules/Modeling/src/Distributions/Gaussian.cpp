#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(unsigned int const dim, double const cov_prec, Gaussian::Mode const mode) : Distribution(), mode(mode), dim(dim), cov(mode==Gaussian::Mode::Covariance? cov_prec : 1.0/cov_prec), prec(mode==Gaussian::Mode::Covariance? 1.0/cov_prec : cov_prec) {}

Gaussian::Gaussian(boost::any const& obj, Gaussian::Mode const mode) : Distribution(), mode(mode), dim(algebra->Size(obj, 0)), cov(mode==Gaussian::Mode::Covariance? SaveCovPrec(obj) : boost::none), prec(mode==Gaussian::Mode::Precision? SaveCovPrec(obj) : boost::none) {}

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
  switch( mode ) {
  case Gaussian::Mode::Covariance:
    return -0.5*algebra->InnerProduct(inputs[0].get(), algebra->ApplyInverse(cov, inputs[0].get()));
  case Gaussian::Mode::Precision:
    return -0.5*algebra->InnerProduct(inputs[0].get(), algebra->Apply(prec, inputs[0].get()));
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
    if( !covSqrt ) {
      covSqrt = algebra->SquareRoot(cov);
    }
    
    // if one dimensional, return a double
    if( dim==1 ) {
      assert(covSqrt->type()==typeid(double));
      return stdnrm(0)*boost::any_cast<double>(*covSqrt);
    }
    
    // otherwise return an Eigen::VectorXd
    return algebra->Apply(*covSqrt, stdnrm);
  }
  case Gaussian::Mode::Precision: {
    if( !precSqrt ) {
      precSqrt = algebra->SquareRoot(prec);
    }

    // if one dimensional, return a double
    if( dim==1 ) {
      assert(precSqrt->type()==typeid(double));
      return stdnrm(0)/boost::any_cast<double>(*precSqrt);
    }
    
    // otherwise return an Eigen::VectorXd
    return algebra->ApplyInverse(*precSqrt, stdnrm);
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
