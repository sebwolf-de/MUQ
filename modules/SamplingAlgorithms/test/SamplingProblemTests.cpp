#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

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
