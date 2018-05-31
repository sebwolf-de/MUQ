/***
# Example 2: Introduction to block MCMC sampling
## Overview
The goal of this example is to demonstrate MCMC sampling with block updates.
The problem is to sample a Gaussian distribution where the variance of
the Gaussian is a random variable endowed with an inverse Gamma distribution.
Thus, there are two "blocks" of interest: the Gaussian random variable and the
variance hyperparameter.  We will sample the joint distribution of both parameters
by constructing a Metropolis-Within-Gibbs sampler.  In this
setting, the Metropolis-Hastings rule is used in one block while the
other block is fixed.

### Problem Formulation
Let $x$ denote the Gaussian random variable of interest, and $\sigma^2$ its prior
variance.  Notice that the joint density $\pi(x,\sigma)$ can be expanded in two
ways:
$$
\pi(x,\sigma^2) = \pi(x | \sigma^2)\pi(\sigma^2)
$$
and
$$
\pi(x,\sigma^2) = \pi(\sigma^2 | x)\pi(x).
$$
This will be useful below.


Now consider a two stage Metropolis-Hastings algorithm.  In the first stage,
we take a step in $x$ for a fixed value of $\sigma^2$.  In the second stage,
we take a step in $\sigma^2$ for a fixed value of $x$.  The algorithm would look
something like this:

1. Update $x$

  a. Propose a step in the $x$ block, $x^\prime \sim q(x | x^{(k-1)}, \sigma^{(k-1)})$

  b. Compute the acceptance probability using the expanded joint density
  $$\begin{eqnarray}
  \gamma &=& \frac{\pi(x^\prime | \sigma^{(k-1)})\pi(\sigma^{(k-1)})}{\pi(x^{(k-1)} | \sigma^{(k-1)}) \pi(\sigma^{(k-1)})} \frac{q(x^{(k-1)} | x^\prime, \sigma^{(k-1)})}{q(x^\prime | x^{(k-1)}, \sigma^{(k-1)})} \\\\
         &=& \frac{\pi(x^\prime | \sigma^{(k-1)})}{\pi(x^{(k-1)} | \sigma^{(k-1)})} \frac{q(x^{(k-1)} | x^\prime, \sigma^{(k-1)})}{q(x^\prime | x^{(k-1)}, \sigma^{(k-1)})}
  \end{eqnarray}$$

  c. Take the step in the $x$ block: $x^{(k)} = x^\prime$ with probability $\min(1,\gamma)$, else $x^{(k)} = x^{(k-1)}$

2. Update $\sigma^2$

  a. Propose a step in the $\sigma^2$ block, $\sigma^\prime \sim q(\sigma | x^{(k)}, \sigma^{(k-1)})$

  b. Compute the acceptance probability using the expanded joint density
  $$\begin{eqnarray}
  \gamma &=& \frac{\pi(\sigma^\prime | x^{(k)})\pi(x^{(k)})}{\pi(\sigma^{(k-1)} | x^{(k)}) \pi(x^{(k)})} \frac{q(\sigma^{(k-1)} | \sigma^\prime, x^{(k)})}{q(\sigma^\prime | \sigma^{(k-1)}, x^{(k)})}. \\\\
         &=& \frac{\pi(\sigma^\prime | x^{(k)})}{\pi(\sigma^{(k-1)} | x^{(k)})} \frac{q(\sigma^{(k-1)} | \sigma^\prime, x^{(k)})}{q(\sigma^\prime | \sigma^{(k-1)}, x^{(k)})}.
  \end{eqnarray}$$

  c. Take the step in the $\sigma^2$ block: $\sigma^{(k)} = \sigma^\prime$ with probability $\min(1,\gamma)$, else $\sigma^{(k)} = \sigma^{(k-1)}$

The extra complexity of this two stage approach is warranted when one or both of the block proposals
$q(\sigma^\prime | \sigma^{(k-1)}, x^{(k)})$ and $q(x^\prime | x^{(k-1)}, \sigma^{(k-1)})$
can be chosen to match the condtiional target densities $\pi(\sigma^\prime | x^{(k)})$
and $\pi(x^\prime | \sigma^{(k-1)})$.  For example, when $\pi(x | \sigma^2)$ is Gaussian
and $\pi(\sigma^2)$ is Inverse Gamma, then $\pi(\sigma^2 | x)$ can be sampled
directly, allowing us to choose $q(\sigma^\prime | \sigma^{(k-1)}, x^{(k)}) = \pi(\sigma^\prime | x^{(k)})$,
guaranteeing an acceptance probability of one for the $\sigma^2$ update.  Notice
that in this example, $\pi(x^\prime | \sigma^{(k-1)})$ is Gaussian and can also
be sampled directly.  However, for illustrative purposes, we will mix a random
walk proposal on $x$ with an independent Inverse Gamma proposal on $\sigma^2$.

*/
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/Modeling/IdentityOperator.h"
#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;


int main(){

  auto varPiece = std::make_shared<IdentityOperator>(1);

  Eigen::VectorXd mu(2);
  mu << 1.0, 2.0;

  auto gaussDens = std::make_shared<Gaussian>(mu, Gaussian::DiagCovariance)->AsDensity();

  std::cout << "Gaussian piece has " << gaussDens->inputSizes.size()
            << " inputs with sizes " << gaussDens->inputSizes.transpose() << std::endl;

  const double alpha = 2.5;
  const double beta = 1.0;


  auto varDens = std::make_shared<InverseGamma>(alpha,beta)->AsDensity();
  auto prodDens = std::make_shared<DensityProduct>(2);

  auto replOp = std::make_shared<ReplicateOperator>(1,2);

  auto graph = std::make_shared<WorkGraph>();

  graph->AddNode(gaussDens, "Gaussian Density");
  graph->AddNode(varPiece, "Variance");
  graph->AddNode(varDens, "Variance Density");
  graph->AddNode(prodDens, "Joint Density");
  graph->AddNode(replOp, "Replicated Variance");

  graph->AddEdge("Variance", 0, "Replicated Variance", 0);
  graph->AddEdge("Replicated Variance", 0, "Gaussian Density", 1);
  graph->AddEdge("Variance", 0, "Variance Density", 0);

  graph->AddEdge("Gaussian Density", 0, "Joint Density", 0);
  graph->AddEdge("Variance Density", 0, "Joint Density", 1);

  graph->Visualize("DensityGraph.pdf");


  auto jointDens = graph->CreateModPiece("Joint Density");

  auto problem = std::make_shared<SamplingProblem>(jointDens);

  pt::ptree pt;
  pt.put("NumSamples", 1e5); // number of MCMC steps
  pt.put("BurnIn", 1e4);
  pt.put("KernelList", "Kernel1,Kernel2"); // Name of block that defines the transition kernel

  pt.put("Kernel1.Method","MHKernel");  // Name of the transition kernel class
  pt.put("Kernel1.Proposal", "MyProposal"); // Name of block defining the proposal distribution
  pt.put("Kernel1.MyProposal.Method", "MHProposal"); // Name of proposal class
  pt.put("Kernel1.MyProposal.ProposalVariance", 0.5); // Variance of the isotropic MH proposal

  pt.put("Kernel2.Method","MHKernel");  // Name of the transition kernel class
  pt.put("Kernel2.Proposal", "GammaProposal"); // Name of block defining the proposal distribution

  pt.put("Kernel2.GaussProposal.Method", "MHProposal");
  pt.put("Kernel2.GaussProposal.ProposalVariance", 0.25);

  pt.put("Kernel2.GammaProposal.Method", "InverseGammaProposal");
  pt.put("Kernel2.GammaProposal.InverseGammaNode", "Variance Density");
  pt.put("Kernel2.GammaProposal.GaussianNode", "Gaussian Density");

  /***
  Once the algorithm parameters are specified, we can pass them to the SingleChainMCMC
  constructor to create an instance of the MCMC algorithm we defined in the
  property tree.
  */
  auto mcmc = std::make_shared<SingleChainMCMC>(pt,problem);


  /***
  ### 3. Run the MCMC algorithm
  We are now ready to run the MCMC algorithm.  Here we start the chain at the
  target densities mean.   The resulting samples are returned in an instance
  of the SampleCollection class, which internally holds the steps in the MCMC chain
  as a vector of weighted SamplingState's.
  */
  std::vector<Eigen::VectorXd> startPt(2);
  startPt.at(0) = mu;
  startPt.at(1) = Eigen::VectorXd::Ones(1);

  SampleCollection const& samps = mcmc->Run(startPt);

  Eigen::VectorXd sampMean = samps.Mean();
  std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;

  Eigen::VectorXd sampVar = samps.Variance();
  std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;

  Eigen::MatrixXd sampCov = samps.Covariance();
  std::cout << "\nSample Covariance = \n" << sampCov << std::endl;

  Eigen::VectorXd sampMom3 = samps.CentralMoment(3);
  std::cout << "\nSample Third Moment = \n" << sampMom3 << std::endl << std::endl;

  return 0;
}
