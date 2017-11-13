#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

TEST(ImportanceSampling, Setup) {
  // create an instance of importance sampling
  auto ip = std::make_shared<ImportanceSampling>();

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", 100); // number of Monte Carlo samples

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>();

  // evaluate
  ip->Evaluate(pt, problem);
}
