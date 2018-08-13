#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/GMHKernel.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(GMHKernelTest, PrposalTest) {
  // create the density we wish to sample
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity();

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // use an apative Metropolis proposal
  pt::ptree pt;
  pt.put("MyKernel.MyProposal.ProposalVariance", 1.0);
  pt.put("MyKernel.MyProposal.AdaptSteps", 200);
  pt.put("MyKernel.MyProposal.AdaptStart", 2000);
  pt.put("MyKernel.MyProposal.AdaptScale", 1.0);
  pt.put("MyKernel.MyProposal.AdaptScale", 1.0);
  pt.put("MyKernel.MyProposal.Method", "AMProposal");
  pt.put("MyKernel.Proposal", "MyProposal");
  pt.put("MyKernel.NumProposals", 10);

  // create the kernel
  auto kern = std::make_shared<GMHKernel>(pt.get_child("MyKernel"), problem);

  // some dummy state
  auto state = std::make_shared<SamplingState>(Eigen::VectorXd::Random(2), 1.0);

  // propose a bunch of points in the pre step
  kern->PreStep(0, state);
}
