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
  pt.put<std::string>("SamplingAlgorithm.TransitionKernel", "MHKernel"); // the transition kernel

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>();

  // evaluate
  mcmc->Evaluate(pt, problem);
}
