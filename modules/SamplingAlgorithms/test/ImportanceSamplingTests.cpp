#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(ImportanceSampling, Setup) {
  // create an instance of importance sampling
  auto ip = std::make_shared<ImportanceSampling>();

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", 100); // number of Monte Carlo samples

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  auto dist = std::make_shared<Gaussian>(); // it is standard normal (1D) by default

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // evaluate
  ip->Evaluate(pt, problem);
}
