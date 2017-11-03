#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

using namespace muq::Modeling;

TEST(GaussianDistributionTests, IsotropicDensity) {
  // create the distributions
  auto standard1D = std::make_shared<Gaussian>();
  auto standard = std::make_shared<Gaussian>(3);

  // compute the log density for a standard normal
  const double x = 2.5;
  const double logstandard1D = standard1D->LogDensity(x);

  const Eigen::Vector3d x3(2.0, 1.4, 3.0);
  const double logstandard = standard->LogDensity(x3);

  EXPECT_DOUBLE_EQ(logstandard1D, -x*x/2.0);
  EXPECT_DOUBLE_EQ(logstandard, -x3.dot(x3)/2.0);

  const unsigned int N = 1.0e6;

  double mean = boost::any_cast<double>(standard1D->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    mean += boost::any_cast<double>(standard1D->Sample());
  }
  mean /= (N+1.0);
  
  EXPECT_NEAR(mean, 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, ScaledIsotropicDensity) {
  // the covariance or precision scale (variance in the cov. case and 1/variance in the prec. case)
  const double sigma2 = 0.25; 

  auto scaledIdentityCov1D = std::make_shared<Gaussian>(1, sigma2);
  auto scaledIdentityCov = std::make_shared<Gaussian>(3, sigma2);
  auto scaledIdentityPrec = std::make_shared<Gaussian>(3, sigma2, Gaussian::Mode::Precision);
    
  const Eigen::Vector3d x(2.0, 1.4, 3.0);
  const double logIdentityCov = scaledIdentityCov->LogDensity(x);
  const double logIdentityPrec = scaledIdentityPrec->LogDensity(x);

  EXPECT_DOUBLE_EQ(logIdentityCov, -x.dot(x)/(2.0*sigma2));
  EXPECT_DOUBLE_EQ(logIdentityPrec, -x.dot(x)*sigma2/2.0);

  const unsigned int N = 1.0e6;
  Eigen::VectorXd samps(N+1);
  
  samps(0) = boost::any_cast<double>(scaledIdentityCov1D->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    samps(i+1) = boost::any_cast<double>(scaledIdentityCov1D->Sample());
  }

  const double mean = samps.sum()/(N+1.0);
  samps = samps.array()-mean;
  const double cov = (samps.array()*samps.array()).sum()/N;

  EXPECT_NEAR(cov, sigma2, 1.0e-2);
  EXPECT_NEAR(mean, 0.0, 1.0e-2);
}
