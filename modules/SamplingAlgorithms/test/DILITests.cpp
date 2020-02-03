#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/LinearAlgebra/HessianOperator.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/DILIKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"

#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"
#include "MUQ/Approximation/GaussianProcesses/MaternKernel.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

std::shared_ptr<Gaussian> DILITest_CreatePrior(unsigned int N)
{
  Eigen::MatrixXd xs = Eigen::VectorXd::LinSpaced(N, 0,1).transpose();

  MaternKernel kern(1, 1.0, 0.05, 5.0/2.0);
  ZeroMean mu(1,1);
  GaussianProcess gpPrior(mu,kern);

  return gpPrior.Discretize(xs);
}

TEST(MCMC, DILIKernel_HessianOperator) {

  const unsigned int numNodes = 100;
  const unsigned int dataDim = 10;
  const double noiseStd = 1e-3;

  std::shared_ptr<Gaussian> prior = DILITest_CreatePrior(numNodes);

  Eigen::MatrixXd forwardMat = Eigen::MatrixXd::Zero(dataDim,numNodes);
  forwardMat.row(0) = Eigen::RowVectorXd::Ones(numNodes);
  for(unsigned int i=1; i<dataDim; ++i)
    forwardMat(i,10*i) = 1.0;

  Eigen::MatrixXd trueHess = (1.0/(noiseStd*noiseStd))*forwardMat.transpose() * forwardMat;

  auto forwardMod = LinearOperator::Create(forwardMat);
  Eigen::VectorXd data = forwardMod->Apply(prior->Sample());

  auto noiseModel = std::make_shared<Gaussian>(data, noiseStd*noiseStd*Eigen::VectorXd::Ones(dataDim));

  WorkGraph graph;
  graph.AddNode(std::make_shared<IdentityOperator>(numNodes), "Parameters");
  graph.AddNode(prior->AsDensity(), "Prior");
  graph.AddNode(forwardMod, "ForwardModel");
  graph.AddNode(noiseModel->AsDensity(), "Likelihood");
  graph.AddNode(std::make_shared<DensityProduct>(2),"Posterior");
  graph.AddEdge("Parameters",0,"Prior",0);
  graph.AddEdge("Parameters",0,"ForwardModel",0);
  graph.AddEdge("ForwardModel",0,"Likelihood",0);
  graph.AddEdge("Prior",0,"Posterior",0);
  graph.AddEdge("Likelihood",0,"Posterior",1);

  auto logLikely = graph.CreateModPiece("Likelihood");

  std::vector<Eigen::VectorXd> inputs{prior->GetMean()};
  Eigen::VectorXd sens = Eigen::VectorXd::Ones(1);
  auto hessOp = std::make_shared<HessianOperator>(logLikely, inputs, 0, 0, 0, sens, -1.0);

  Eigen::MatrixXd hess(numNodes,numNodes);
  for(unsigned int i=0; i<numNodes; ++i){
    Eigen::VectorXd unitVec = Eigen::VectorXd::Zero(numNodes);
    unitVec(i) = 1.0;
    hess.col(i) = hessOp->Apply(unitVec);
  }

  for(unsigned int j=0;j<numNodes;++j){
    for(unsigned int i=0; i<numNodes; ++i){
      EXPECT_NEAR(trueHess(i,j), hess(i,j), 1e-5*std::abs(trueHess(i,j)));
    }
  }
}


TEST(MCMC, DILIKernel_Construction) {

  const unsigned int numNodes = 100;
  const unsigned int dataDim = 5;
  const double noiseStd = 1e-2;

  std::shared_ptr<Gaussian> prior = DILITest_CreatePrior(numNodes);

  Eigen::MatrixXd forwardMat = Eigen::MatrixXd::Zero(dataDim,numNodes);
  forwardMat.row(0) = Eigen::RowVectorXd::Ones(numNodes)/numNodes;
  for(unsigned int i=1; i<dataDim; ++i)
    forwardMat(i,10*i) = 1.0;

  auto forwardMod = LinearOperator::Create(forwardMat);
  Eigen::VectorXd trueField = Eigen::VectorXd::Zero(numNodes);
  Eigen::VectorXd data = forwardMod->Apply(trueField);

  auto noiseModel = std::make_shared<Gaussian>(data, noiseStd*noiseStd*Eigen::VectorXd::Ones(dataDim));

  WorkGraph graph;
  graph.AddNode(std::make_shared<IdentityOperator>(numNodes), "Parameters");
  graph.AddNode(prior->AsDensity(), "Prior");
  graph.AddNode(forwardMod, "ForwardModel");
  graph.AddNode(noiseModel->AsDensity(), "Likelihood");
  graph.AddNode(std::make_shared<DensityProduct>(2),"Posterior");
  graph.AddEdge("Parameters",0,"Prior",0);
  graph.AddEdge("Parameters",0,"ForwardModel",0);
  graph.AddEdge("ForwardModel",0,"Likelihood",0);
  graph.AddEdge("Prior",0,"Posterior",0);
  graph.AddEdge("Likelihood",0,"Posterior",1);

  auto logLikely = graph.CreateModPiece("Likelihood");

  pt::ptree pt;
  const unsigned int numSamps = 10000;
  pt.put("NumSamples",numSamps);
  pt.put("BurnIn", 0);
  pt.put("PrintLevel",3);
  pt.put("HessianType","Exact");

  pt.put("Eigensolver Block", "LOBPCG");
  pt.put("LOBPCG.NumEigs",15); // Maximum number of generalized eigenvalues to compute (e.g., maximum LIS dimension)
  pt.put("LOBPCG.RelativeTolerance", 1e-3); // Fraction of the largest eigenvalue used as stopping criteria on how many eigenvalues to compute
  pt.put("LOBPCG.AbsoluteTolerance",0.0); // Minimum allowed eigenvalue
  pt.put("LOBPCG.BlockSize",20);
  pt.put("LOBPCG.MaxIts",200);
  pt.put("LOBPCG.Verbosity",0);
  pt.put("LOBPCG.SolverTolerance",1e-6);

  pt.put("LIS Block", "LIS");
  pt.put("LIS.Method", "MHKernel");
  pt.put("LIS.Proposal","MyProposal");
  pt.put("LIS.MyProposal.Method","MHProposal");
  pt.put("LIS.MyProposal.ProposalVariance", 1.0);

  pt.put("CS Block", "CS");
  pt.put("CS.Method", "MHKernel");
  pt.put("CS.Proposal","MyProposal");
  pt.put("CS.MyProposal.Method", "CrankNicolsonProposal");
  pt.put("CS.MyProposal.Beta",1.0);
  pt.put("CS.MyProposal.PriorNode","Prior");

  // The posterior is Gaussian in this case, so we can analytically compute the posterior mean and covariance for comparison
  std::shared_ptr<Gaussian> truePost = prior->Condition(forwardMat, data, noiseStd*noiseStd*Eigen::VectorXd::Ones(dataDim));

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);


  auto kernel = std::make_shared<DILIKernel>(pt, problem, logLikely, prior);

  auto sampler = std::make_shared<SingleChainMCMC>(pt, std::vector<std::shared_ptr<TransitionKernel>>{kernel});

  auto samps = sampler->Run(truePost->Sample()); // Use a true posterior sample to avoid burnin


  // Check the posterior moments
  Eigen::VectorXd sampMean = samps->Mean();
  Eigen::MatrixXd sampCov = samps->Covariance();
  Eigen::VectorXd trueMean = truePost->GetMean();
  Eigen::MatrixXd trueCov = truePost->GetCovariance();

  Eigen::VectorXd ess = samps->ESS();

  for(int i=0; i<numNodes; ++i){

    // Estimator variance
    double estVar = trueCov(i,i)/ess(i);
    EXPECT_NEAR(trueMean(i), sampMean(i), 3.0*std::sqrt(estVar));
  }
}
