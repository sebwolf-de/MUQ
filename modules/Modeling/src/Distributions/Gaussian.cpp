#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(unsigned int const dim, double const cov_prec, Gaussian::Mode const mode) : Distribution(), mode(mode), dim(dim), cov(mode==Gaussian::Mode::Covariance? cov_prec : 1.0/cov_prec), prec(mode==Gaussian::Mode::Covariance? 1.0/cov_prec : cov_prec) {}

Gaussian::Gaussian(boost::any const& diag, Gaussian::Mode const mode) : Distribution(), mode(mode), dim(algebra->Size(diag)), cov(mode==Gaussian::Mode::Covariance? diag : boost::none), prec(mode==Gaussian::Mode::Precision? diag : boost::none) {}

Gaussian::~Gaussian() {}

double Gaussian::LogDensityImpl(ref_vector<boost::any> const& inputs) const {
  switch( mode ) {
  case Gaussian::Mode::Covariance: 
    return -0.5*algebra->InnerProduct(inputs[0].get(), algebra->ApplyInverse(cov, inputs[0].get()));
  case Gaussian::Mode::Precision:
    return -0.5*algebra->InnerProduct(inputs[0].get(), algebra->Apply(prec, inputs[0].get()));
  default:
    assert(false);
  }
}

boost::any Gaussian::SampleImpl(ref_vector<boost::any> const& inputs) const {
  // make sure the dimension is nonzero
  assert(dim>0);

  const Eigen::VectorXd stdnrm = RandomGenerator::GetNormal(dim);
  
  // if one dimensional, return a double
  if( dim==1 ) {
    return stdnrm(0)*std::sqrt(boost::any_cast<double>(cov));
  }

  // otherwise return an Eigen::VectorXd
  return stdnrm*std::sqrt(boost::any_cast<double>(cov));
}
