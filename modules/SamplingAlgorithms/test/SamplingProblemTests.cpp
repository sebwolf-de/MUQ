#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/UniformBox.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(SamplingProblem, GaussianTarget) {
  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  auto dist = std::make_shared<Gaussian>(); // it is standard normal (1D) by default

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  const unsigned int N = 1.0e5;
  double mean = boost::any_cast<double const>(problem->SampleTarget(ref_vector<boost::any>()));
  for( unsigned int i=0; i<N; ++i ) {
    mean += boost::any_cast<double const>(problem->SampleTarget(ref_vector<boost::any>()));
  }
  mean /= (double)N;
  
  EXPECT_NEAR(mean, 0.0, 1.0e-2);
}

TEST(SamplingProblem, ImportanceSamplingSetup) {
  // create a Gaussian distribution---this is the biasing distribution
  auto bias = std::make_shared<Gaussian>(); // it is standard normal (1D) by default

  // create a uniform distribution---the sampling problem is built around characterizing this distribution
  auto dist = std::make_shared<UniformBox>(std::pair<double, double>(-0.5, 1.0));

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist, bias);

  // make sure we are correctly sampling the baising distribution
  const unsigned int N = 1.0e5;
  double mean = boost::any_cast<double const>(problem->SampleBiasingDistribution(ref_vector<boost::any>()));
  for( unsigned int i=0; i<N; ++i ) {
    mean += boost::any_cast<double const>(problem->SampleBiasingDistribution(ref_vector<boost::any>()));
  }
  mean /= (double)N;
  EXPECT_NEAR(mean, 0.0, 1.0e-2);

  // make sure we are correctly evaluating the target
  boost::any x = 0.0;
  double logTarget = problem->EvaluateLogTarget(ref_vector<boost::any>(1, std::cref(x)));
  double logBias = problem->EvaluateLogBiasingDistribution(ref_vector<boost::any>(1, std::cref(x)));
  EXPECT_DOUBLE_EQ(logTarget, 1.0);
  EXPECT_DOUBLE_EQ(logBias, bias->LogDensity(x));
  //x = 2.0;
  x = 1.27591;
  logTarget = problem->EvaluateLogTarget(ref_vector<boost::any>(1, std::cref(x)));
  logBias = problem->EvaluateLogBiasingDistribution(ref_vector<boost::any>(1, std::cref(x)));
  EXPECT_TRUE(std::isinf(logTarget));
  EXPECT_DOUBLE_EQ(logBias, bias->LogDensity(x));
}
