#include <gtest/gtest.h>

#include <Eigen/Core>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

TEST(GaussianDistributionTests, SpecifyMean) {
  // create the distributions
  auto standard1D = std::make_shared<Gaussian>();
  const Eigen::Vector3d mu = Eigen::Vector3d::Ones();
  auto standard = std::make_shared<Gaussian>(mu);
  EXPECT_EQ(standard1D->Dimension(), 1);
  EXPECT_EQ(standard->Dimension(), 3);

  // compute the log density for a standard normal
  const double x = 2.5;
  const Eigen::Vector3d x3(2.0, 1.4, 3.0);
  EXPECT_DOUBLE_EQ(standard1D->LogDensity(x), -x*x/2.0);
  EXPECT_DOUBLE_EQ(standard->LogDensity(x3), -(x3-mu).dot(x3-mu)/2.0);

  const unsigned int N = 1.0e5;

  double mean1D = boost::any_cast<double>(standard1D->Sample());
  Eigen::VectorXd mean = boost::any_cast<Eigen::VectorXd>(standard->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    mean1D += boost::any_cast<double>(standard1D->Sample());
    mean += boost::any_cast<Eigen::VectorXd>(standard->Sample());
  }
  mean1D /= (N+1.0);
  mean /= (N+1.0);

  EXPECT_NEAR(mean1D, 0.0, 1.0e-2);
  EXPECT_NEAR((mean-mu).norm(), 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, SpecifyBoth) {
  // create the distributions
  const Eigen::Vector3d cov = Eigen::Vector3d::Random().cwiseAbs();
  const Eigen::Vector3d prec = Eigen::Vector3d::Random().cwiseAbs();
  const Eigen::Vector3d mu = Eigen::Vector3d::Ones();
  auto covDist = std::make_shared<Gaussian>(mu, cov);
  auto precDist = std::make_shared<Gaussian>(mu, prec, Gaussian::Mode::Precision);
  EXPECT_EQ(covDist->Dimension(), 3);
  EXPECT_EQ(precDist->Dimension(), 3);

  // compute the log density for a standard normal
  const Eigen::Vector3d x(2.0, 1.4, 3.0);
  EXPECT_DOUBLE_EQ(covDist->LogDensity(x), -(x-mu).dot((1.0/cov.array()).matrix().asDiagonal()*(x-mu))/2.0);
  EXPECT_DOUBLE_EQ(precDist->LogDensity(x), -(x-mu).dot(prec.asDiagonal()*(x-mu))/2.0);

  const unsigned int N = 2.0e5;

  Eigen::VectorXd meanCov = boost::any_cast<Eigen::VectorXd>(covDist->Sample());
  Eigen::VectorXd meanPrec = boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    meanCov += boost::any_cast<Eigen::VectorXd>(covDist->Sample());
    meanPrec += boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  }
  meanCov /= (N+1.0);
  meanPrec /= (N+1.0);

  EXPECT_NEAR((meanCov-mu).norm(), 0.0, 1.0e-2);
  EXPECT_NEAR((meanPrec-mu).norm(), 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, Scalar) {
  // the covariance or precision scale (variance in the cov. case and 1/variance in the prec. case)
  const double sigma2 = 0.25; 

  auto scalarCov = std::make_shared<Gaussian>(sigma2, Gaussian::Mode::Covariance);
  auto scalarPrec = std::make_shared<Gaussian>(sigma2, Gaussian::Mode::Precision);
  EXPECT_EQ(scalarCov->Dimension(), 1);
  EXPECT_EQ(scalarPrec->Dimension(), 1);

  // test log density
  const double x = 2.25;
  EXPECT_DOUBLE_EQ(scalarCov->LogDensity(x), -x*x/(2.0*sigma2));
  EXPECT_DOUBLE_EQ(scalarPrec->LogDensity(x), -x*x*sigma2/2.0);

  const unsigned int N = 1.0e5;
  Eigen::VectorXd sampsCov(N+1);
  Eigen::VectorXd sampsPrec(N+1);
  
  sampsCov(0) = boost::any_cast<double>(scalarCov->Sample());
  sampsPrec(0) = boost::any_cast<double>(scalarPrec->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    sampsCov(i+1) = boost::any_cast<double>(scalarCov->Sample());
    sampsPrec(i+1) = boost::any_cast<double>(scalarPrec->Sample());
  }

  const double meanCov = sampsCov.sum()/(N+1.0);
  const double meanPrec = sampsPrec.sum()/(N+1.0);
  sampsCov = sampsCov.array()-meanCov;
  sampsPrec = sampsPrec.array()-meanPrec;
  const double cov0 = (sampsCov.array()*sampsCov.array()).sum()/N;
  const double cov1 = N/(sampsPrec.array()*sampsPrec.array()).sum();

  EXPECT_NEAR(cov0, sigma2, 1.0e-2);
  EXPECT_NEAR(cov1, sigma2, 1.0e-2);
  EXPECT_NEAR(meanCov, 0.0, 1.0e-2);
  EXPECT_NEAR(meanPrec, 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, DiagonalCovPrec) {
  // the covariance or precision scale (variance in the cov. case and 1/variance in the prec. case)
  const Eigen::Vector2d covDiag = 1e-2 * Eigen::Vector2d::Ones() + RandomGenerator::GetUniform(2);
  const Eigen::Vector2d precDiag = 1e-2 * Eigen::Vector2d::Ones() + RandomGenerator::GetUniform(2);
  
  auto diagCov = std::make_shared<Gaussian>(covDiag, Gaussian::Mode::Covariance);
  auto diagPrec = std::make_shared<Gaussian>(precDiag, Gaussian::Mode::Precision);
  EXPECT_EQ(diagCov->Dimension(), 2);
  EXPECT_EQ(diagPrec->Dimension(), 2);

  const Eigen::Vector2d x = Eigen::Vector2d::Random();
  EXPECT_DOUBLE_EQ(diagCov->LogDensity(x), -x.dot((1.0/covDiag.array()).matrix().asDiagonal()*x)/2.0);
  EXPECT_DOUBLE_EQ(diagPrec->LogDensity(x), -x.dot(precDiag.asDiagonal()*x)/2.0);
    
  const unsigned int N = 5.0e5;
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
  cov = Eigen::MatrixXd::Identity(dim, dim) + cov*cov.transpose();
  Eigen::MatrixXd prec = Eigen::MatrixXd::Random(dim, dim);
  prec = Eigen::MatrixXd::Identity(dim, dim) + prec*prec.transpose();
  
  Eigen::LLT<Eigen::MatrixXd> covChol;
  covChol.compute(cov);
  const Eigen::MatrixXd covL = covChol.matrixL();

  auto covDist = std::make_shared<Gaussian>(cov, Gaussian::Mode::Covariance);
  auto precDist = std::make_shared<Gaussian>(prec, Gaussian::Mode::Precision);
  EXPECT_EQ(covDist->Dimension(), dim);

  const Eigen::VectorXd x = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd delta = covL.triangularView<Eigen::Lower>().solve(x);
  covL.triangularView<Eigen::Lower>().transpose().solveInPlace(delta);
  EXPECT_DOUBLE_EQ(covDist->LogDensity(x), -x.dot(delta)/2.0);
  EXPECT_DOUBLE_EQ(precDist->LogDensity(x), -x.dot(prec*x)/2.0);

  const unsigned int N = 5.0e5;
  Eigen::VectorXd meanCov = boost::any_cast<Eigen::VectorXd>(covDist->Sample());
  Eigen::VectorXd meanPrec = boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    meanCov += boost::any_cast<Eigen::VectorXd>(covDist->Sample());
    meanPrec += boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  }
  meanCov /= (N+1.0);
  meanPrec /= (N+1.0);

  EXPECT_NEAR(meanCov.norm(), 0.0, 1.0e-2);
  EXPECT_NEAR(meanPrec.norm(), 0.0, 1.0e-2);
}

