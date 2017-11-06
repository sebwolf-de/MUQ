#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(unsigned int const dim, double const cov_prec, Gaussian::Mode const mode) : Distribution(), dim(dim), cov(mode==Gaussian::Mode::Covariance? cov_prec : 1.0/cov_prec) {
  // initialize the any algebra
  algebra = std::make_shared<AnyAlgebra>();
}

Gaussian::~Gaussian() {}

double Gaussian::LogDensityImpl(ref_vector<boost::any> const& inputs) const {
  return -0.5*algebra->InnerProduct(inputs[0].get(), inputs[0].get())/cov;
}

boost::any Gaussian::SampleImpl(ref_vector<boost::any> const& inputs) const {
  // make sure the dimension is nonzero
  assert(dim>0);

  const Eigen::VectorXd stdnrm = RandomGenerator::GetNormal(dim);
  
  // if one dimensional, return a double
  if( dim==1 ) {
    return stdnrm(0)*sqrt(cov);
  }

  // otherwise return an Eigen::VectorXd
  return stdnrm*sqrt(cov);
}
