#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/DRKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/MALAProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/Utilities/AnyHelpers.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

TEST(MCMC, DRKernel_MHProposal) {
  const unsigned int N = 1e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("PrintLevel",0);
  pt.put("KernelList", "Kernel1"); // the transition kernel
  pt.put("Kernel1.Method","DRKernel");
  pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel1.MyProposal.Method", "MHProposal");
  pt.put("Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal
  pt.put("Kernel1.NumStages", 3);
  pt.put("Kernel1.ScaleFunction", "Linear"); // the variance of the isotropic MH proposal
  pt.put("Kernel1.Scale", 1.0);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt,problem);

  // Make sure the kernel and proposal are correct
  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  ASSERT_TRUE(kernelBase);
  std::shared_ptr<DRKernel> kernelDR = std::dynamic_pointer_cast<DRKernel>(kernelBase);
  ASSERT_TRUE(kernelDR);

  std::vector<std::shared_ptr<MCMCProposal>> proposals = kernelDR->Proposals();
  EXPECT_EQ(pt.get<int>("Kernel1.NumStages"), proposals.size());

  std::shared_ptr<MHProposal> proposalMH = std::dynamic_pointer_cast<MHProposal>(proposals.at(0));
  ASSERT_TRUE(proposalMH);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  EXPECT_EQ(pt.get<int>("NumSamples"), samps->size());

  //boost::any anyMean = samps.Mean();
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}
