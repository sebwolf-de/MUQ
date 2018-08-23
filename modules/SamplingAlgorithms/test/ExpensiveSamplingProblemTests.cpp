#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(ExpensiveSamplingProblemTests, GaussianTarget) {
  const unsigned int N = 5e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  pt.put("MySamplingProblem.RegressionOptions", "MyRegression");
  pt.put("MySamplingProblem.MyRegression.NumNeighbors", 5);
  pt.put("MySamplingProblem.MyRegression.Order", 1); // approximating the quardatic log-Gaussian with a locally linear function

  pt.put("MySamplingProblem.BetaScale", 1.0);
  pt.put("MySamplingProblem.BetaExponent", 0.9);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<ExpensiveSamplingProblem>(dist, pt.get_child("MySamplingProblem"));

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem);

  // run MCMC
  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  Eigen::VectorXd mean = samps->Mean();

  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}
