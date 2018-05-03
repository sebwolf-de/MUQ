#include <gtest/gtest.h>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/AnyHelpers.h"

#include "MUQ/SamplingAlgorithms/SampleCollection.h"


using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

class SampleCollectionTest : public::testing::Test {
protected:

    virtual void SetUp() override {

      L.resize(2,2);
      L << 1.0, 0.0,
          1.0, 2.0;

      samps = L * RandomGenerator::GetNormal(2,numSamps);
      weights = RandomGenerator::GetUniform(numSamps);

      for(int i=0; i<numSamps; ++i)
        collection.Add(std::make_shared<SamplingState>(Eigen::VectorXd(samps.col(i)), 1.0));

    }

    virtual void TearDown() override {
    }

    const int numSamps = 1e4;
    Eigen::MatrixXd L;

    Eigen::MatrixXd samps;
    Eigen::VectorXd weights;
    SampleCollection collection;
};

TEST_F(SampleCollectionTest, Mean)
{
  boost::any muAny = collection.Mean(0);

  Eigen::VectorXd const& mu = AnyConstCast(muAny);
  Eigen::VectorXd trueMu = samps.rowwise().mean();

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 3*mcStd);

  EXPECT_NEAR(trueMu(0), mu(0), 1e-13);
  EXPECT_NEAR(trueMu(1), mu(1), 1e-13);

  muAny = collection.Mean();
  Eigen::VectorXd const& mu2 = AnyConstCast(muAny);
  mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu2(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu2(1), 3*mcStd);

  EXPECT_NEAR(trueMu(0), mu2(0), 1e-13);
  EXPECT_NEAR(trueMu(1), mu2(1), 1e-13);
}

TEST_F(SampleCollectionTest, SampMean)
{
  boost::any muAny = collection.Mean(0);

  Eigen::VectorXd mu = samps.rowwise().mean();

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 3*mcStd);
}


TEST_F(SampleCollectionTest, Variance)
{
  boost::any anyVar = collection.Variance(0);

  Eigen::VectorXd sampMu = samps.rowwise().mean();
  Eigen::VectorXd sampVar = (samps.colwise() - sampMu).array().pow(2.0).matrix().rowwise().mean();

  Eigen::VectorXd const& var = AnyCast(anyVar);
  EXPECT_NEAR(L(0,0)*L(0,0), var(0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(L(1,0)*L(1,0) + L(1,1)*L(1,1), var(1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampVar(0), var(0), 1e-13);
  EXPECT_NEAR(sampVar(1), var(1), 1e-13);

  anyVar = collection.Variance();
  Eigen::VectorXd const& var2 = AnyCast(anyVar);
  EXPECT_NEAR(L(0,0)*L(0,0), var2(0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(L(1,0)*L(1,0) + L(1,1)*L(1,1), var2(1), 50.0/sqrt(double(numSamps)));

}

TEST_F(SampleCollectionTest, Covariance)
{
  boost::any anyCov = collection.Covariance(0);
  Eigen::MatrixXd const& cov = AnyConstCast(anyCov);

  Eigen::VectorXd sampMu  = samps.rowwise().mean();
  Eigen::MatrixXd sampCov = (samps.colwise()-sampMu) * (samps.colwise()-sampMu).transpose() / (samps.cols());

  Eigen::MatrixXd trueCov = L*L.transpose();

  EXPECT_NEAR(trueCov(0,0), cov(0,0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), cov(0,1), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), cov(1,0), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), cov(1,1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampCov(0,0), cov(0,0), 1e-13);
  EXPECT_NEAR(sampCov(0,1), cov(0,1), 1e-13);
  EXPECT_NEAR(sampCov(1,0), cov(1,0), 1e-13);
  EXPECT_NEAR(sampCov(1,1), cov(1,1), 1e-13);

  anyCov = collection.Covariance();
  Eigen::MatrixXd const& cov2 = AnyConstCast(anyCov);

  EXPECT_NEAR(trueCov(0,0), cov2(0,0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), cov2(0,1), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), cov2(1,0), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), cov2(1,1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampCov(0,0), cov2(0,0), 1e-13);
  EXPECT_NEAR(sampCov(0,1), cov2(0,1), 1e-13);
  EXPECT_NEAR(sampCov(1,0), cov2(1,0), 1e-13);
  EXPECT_NEAR(sampCov(1,1), cov2(1,1), 1e-13);
}
