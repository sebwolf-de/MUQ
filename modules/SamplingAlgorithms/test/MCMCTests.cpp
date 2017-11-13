#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MCMC.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

TEST(MCMC, Setup) {
  // create an instance of MCMC
  auto mcmc = std::make_shared<MCMC>();

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", 100); // number of Monte Carlo samples

  // evaluate
  mcmc->Evaluate(pt);
}
