#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/MCMC.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(MCMC, MHKernel_MHProposal) {
  // create an instance of MCMC
  auto mcmc = std::make_shared<MCMC>();

  const unsigned int N = 7.0e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", N); // number of Monte Carlo samples
  pt.put<std::string>("SamplingAlgorithm.TransitionKernel", "MHKernel"); // the transition kernel
  pt.put<std::string>("MCMC.Proposal", "MHProposal"); // the proposal
  pt.put<unsigned int>("MCMC.StateDimension", 2); // the state dimension
  pt.put<double>("MCMC.Proposal.MH.ProposalSize", 0.75); // the size of the MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // evaluate
  mcmc->Evaluate(pt, problem, start);

  // estimate the mean
  const boost::any mean = mcmc->FirstMoment();
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXd const>(mean)-mu).norm(), 0.0, 1.0e-2);
}

TEST(MCMC, MHKernel_AMProposal) {
  // create an instance of MCMC
  auto mcmc = std::make_shared<MCMC>();

  const unsigned int N = 2.0e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", N); // number of Monte Carlo samples
  pt.put<std::string>("SamplingAlgorithm.TransitionKernel", "MHKernel"); // the transition kernel
  pt.put<std::string>("MCMC.Proposal", "AMProposal"); // the proposal
  pt.put<unsigned int>("MCMC.StateDimension", 2); // the state dimension
  pt.put<double>("MCMC.Proposal.MH.ProposalSize", 0.75); // the size of the MH proposal
  pt.put<unsigned int>("MCMC.Proposal.AM.Steps", 200); // the number of steps between each adaptation
  //pt.put<unsigned int>("MCMC.Proposal.AM.Start", 200); // start adapting after this step

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // evaluate
  mcmc->Evaluate(pt, problem, start);

  const boost::any mean = mcmc->FirstMoment();
  // estimate the mean
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXd const>(mean)-mu).norm(), 0.0, 1.0e-2);
}
