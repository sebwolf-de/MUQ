#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include <MUQ/Utilities/HDF5/HDF5File.h>

#include <MUQ/Modeling/Distributions/Density.h>

#include <MUQ/SamplingAlgorithms/SingleChainMCMC.h>

#include "Target.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

/**
   @param[in] N The number of MCMC steps
   @pram[in] filename File where we want to store the result, defaults to "" which does not write to file
   @pram[in] dataset The dataset in the file
   \return The effective sample size
 */
double RunMHMCMC(unsigned int const N, std::string const& filename = "") {  
  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel", 1);
  pt.put("MyMCMC.KernelList", "Kernel"); // the transition kernel
  pt.put("MyMCMC.Kernel.Method","MHKernel");
  pt.put("MyMCMC.Kernel.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal
  
  // the target density
  auto dist = std::make_shared<Target>(1.0)->AsDensity(); 
  
  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);
  
  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Zero(1);
  
  // create mcmc
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem);

  // run mcmc
  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);
  if( !filename.empty() ) { samps->WriteToFile(filename); }
  return samps->ESS() (0);
}

int main(int argc, char **argv) {
  const int ierr = MPI_Init(nullptr, nullptr);
    
  const unsigned int N = 1.0e5;
  const unsigned int Ntests = 100;
  const std::string filename = "mh-mcmc.h5";

  double ess = RunMHMCMC(N, filename);

  MPI_Finalize();
}
