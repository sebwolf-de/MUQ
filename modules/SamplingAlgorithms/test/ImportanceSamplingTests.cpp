#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/UniformBox.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(ImportanceSampling, Setup) {
  // create an instance of importance sampling
  auto ip = std::make_shared<ImportanceSampling>();

  // the number of samples
  const unsigned int N = 1.0e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", N); // number of Monte Carlo samples

  // create a Gaussian distribution---this is the biasing distribution
  auto bias = std::make_shared<Gaussian>(0.25); 

  // create a uniform distribution---the sampling problem is built around characterizing this distribution
  auto dist = std::make_shared<UniformBox>(std::pair<double, double>(-0.5, 1.0));

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist, bias);

  // evaluate
  const std::vector<boost::any>& result = ip->Evaluate(pt, problem);
  const std::vector<std::shared_ptr<SamplingState> >& samples = boost::any_cast<std::vector<std::shared_ptr<SamplingState> > const&>(result[0]);
  EXPECT_EQ(samples.size(), N);

  // estimate the mean
  const boost::any mean = ip->FirstMoment();

  EXPECT_NEAR(boost::any_cast<double const>(mean), 0.25, 1.0e-2);
}
