#include <gtest/gtest.h>

#include <Eigen/Core>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

TEST(GaussianDistributionTests, IsotropicDensity) {
  // create the distributions
  auto standard1D = std::make_shared<Gaussian>();
  auto standard = std::make_shared<Gaussian>(3);
  EXPECT_EQ(standard1D->Dimension(), 1);
  EXPECT_EQ(standard->Dimension(), 3);

  // compute the log density for a standard normal
  const double x = 2.5;
  const double logstandard1D = standard1D->LogDensity(x);

  const Eigen::Vector3d x3(2.0, 1.4, 3.0);
  const double logstandard = standard->LogDensity(x3);

  EXPECT_DOUBLE_EQ(logstandard1D, -x*x/2.0);
  EXPECT_DOUBLE_EQ(logstandard, -x3.dot(x3)/2.0);

  const unsigned int N = 1.0e6;

  double mean1D = boost::any_cast<double>(standard1D->Sample());
  Eigen::VectorXd mean = boost::any_cast<Eigen::VectorXd>(standard->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    mean1D += boost::any_cast<double>(standard1D->Sample());
    mean += boost::any_cast<Eigen::VectorXd>(standard->Sample());
  }
  mean1D /= (N+1.0);
  mean /= (N+1.0);

  EXPECT_NEAR(mean1D, 0.0, 1.0e-2);
  EXPECT_NEAR(mean.norm(), 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, ScaledIsotropicDensity) {
  // the covariance or precision scale (variance in the cov. case and 1/variance in the prec. case)
  const double sigma2 = 0.25; 

  auto scaledIdentityCov1D = std::make_shared<Gaussian>(1, sigma2);
  auto scaledIdentityCov = std::make_shared<Gaussian>(3, sigma2);
  auto scaledIdentityPrec = std::make_shared<Gaussian>(3, sigma2, Gaussian::Mode::Precision);
  EXPECT_EQ(scaledIdentityCov1D->Dimension(), 1);
  EXPECT_EQ(scaledIdentityCov->Dimension(), 3);
  EXPECT_EQ(scaledIdentityPrec->Dimension(), 3);
    
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

TEST(GaussianDistributionTests, DiagonalCovPrec) {
  // the covariance or precision scale (variance in the cov. case and 1/variance in the prec. case)
  const Eigen::Vector2d covDiag = 1e-2 * Eigen::Vector2d::Ones() + RandomGenerator::GetUniform(2);
  const Eigen::Vector2d precDiag = 1e-2 * Eigen::Vector2d::Ones() + RandomGenerator::GetUniform(2);
  
  auto diagCov = std::make_shared<Gaussian>(covDiag);
  auto diagPrec = std::make_shared<Gaussian>(precDiag, Gaussian::Mode::Precision);
  EXPECT_EQ(diagCov->Dimension(), 2);
  EXPECT_EQ(diagPrec->Dimension(), 2);

  const Eigen::Vector2d x = Eigen::Vector2d::Random();
  EXPECT_DOUBLE_EQ(diagCov->LogDensity(x), -x.dot((1.0/covDiag.array()).matrix().asDiagonal()*x)/2.0);
  EXPECT_DOUBLE_EQ(diagPrec->LogDensity(x), -x.dot(precDiag.asDiagonal()*x)/2.0);
    
  const unsigned int N = 1.0e6;
  Eigen::Vector2d meanCov = boost::any_cast<Eigen::VectorXd>(diagCov->Sample());
  Eigen::Vector2d meanPrec = boost::any_cast<Eigen::VectorXd>(diagCov->Sample());
  
  for( unsigned int i=0; i<N; ++i ) {
    meanCov += boost::any_cast<Eigen::VectorXd>(diagCov->Sample());
    meanPrec += boost::any_cast<Eigen::VectorXd>(diagPrec->Sample());
  }
  meanCov /= (N+1.0);
  meanPrec /= (N+1.0);

  EXPECT_NEAR(meanCov.norm(), 0.0, 1.0e-2);
  EXPECT_NEAR(meanPrec.norm(), 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, MatrixCovPrec) {
  const unsigned int dim = 5;
  Eigen::MatrixXd cov = Eigen::MatrixXd::Random(dim, dim);
  cov = 1e-4*Eigen::MatrixXd::Identity(dim, dim) + cov*cov.transpose();
  Eigen::MatrixXd prec = Eigen::MatrixXd::Random(dim, dim);
  prec = 1e-4*Eigen::MatrixXd::Identity(dim, dim) + prec*prec.transpose();
  
  Eigen::LLT<Eigen::MatrixXd> covChol;
  covChol.compute(cov);
  const Eigen::MatrixXd covL = covChol.matrixL();

  Eigen::LLT<Eigen::MatrixXd> precChol;
  precChol.compute(prec);
  const Eigen::MatrixXd precL = precChol.matrixL();

  auto covDist = std::make_shared<Gaussian>(cov);
  auto precDist = std::make_shared<Gaussian>(prec, Gaussian::Mode::Precision);
  EXPECT_EQ(covDist->Dimension(), dim);

  const Eigen::VectorXd x = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd delta = covL.triangularView<Eigen::Lower>().solve(x);
  covL.triangularView<Eigen::Lower>().transpose().solveInPlace(delta);
  EXPECT_DOUBLE_EQ(covDist->LogDensity(x), -x.dot(delta)/2.0);

  //EXPECT_DOUBLE_EQ(precDist->LogDensity(x), -x.dot(prec*x)/2.0);
}
