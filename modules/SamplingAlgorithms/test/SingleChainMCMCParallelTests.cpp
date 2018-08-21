#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

TEST(MCMC, MHKernel_MHProposal) {
  const unsigned int N = 1e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"),problem);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);
  
  auto comm = mcmc->GetCommunicator();

  if( comm->GetRank()==0 ) {
    EXPECT_EQ(pt.get<int>("MyMCMC.NumSamples"), samps->size());
    
    //boost::any anyMean = samps.Mean();
    Eigen::VectorXd mean = samps->Mean();
    EXPECT_NEAR(mu(0), mean(0), 1e-1);
    EXPECT_NEAR(mu(1), mean(1), 1e-1);
    
    Eigen::MatrixXd cov = samps->Covariance();
    EXPECT_NEAR(1.0, cov(0,0), 1e-1);
    EXPECT_NEAR(0.0, cov(0,1), 1e-1);
    EXPECT_NEAR(0.0, cov(1,0), 1e-1);
    EXPECT_NEAR(1.0, cov(1,1), 1e-1);
  } else {
    EXPECT_EQ(0, samps->size());
  }
}