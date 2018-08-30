#include <boost/property_tree/ptree.hpp>

#include <MUQ/Utilities/HDF5/HDF5File.h>

#include <MUQ/Modeling/Distributions/Density.h>

#include <MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h>
#include <MUQ/SamplingAlgorithms/SingleChainMCMC.h>

#include "Target.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

/**
   @param[in] N The number of MCMC steps
   @pram[in] filename File where we want to store the result, defaults to "" which does not write to file
 */
std::vector<Eigen::MatrixXd> RunLAMCMC(unsigned int const N, std::string const& dataset = "", std::string const& filename = "") {  
  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel", 1);
  pt.put("MyMCMC.KernelList", "Kernel"); // the transition kernel
  pt.put("MyMCMC.Kernel.Method","MHKernel");
  pt.put("MyMCMC.Kernel.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  pt.put("MySamplingProblem.RegressionOptions", "MyRegression");
  pt.put("MySamplingProblem.MyRegression.NumNeighbors", 8);
  //pt.put("MySamplingProblem.MyRegression.Order", 1);
  pt.put("MySamplingProblem.MyRegression.Order", 5);

  pt.put("MySamplingProblem.StructuralScaling", 1.0);
  pt.put("MySamplingProblem.PoisednessConstant", 1.0);
  //pt.put("MySamplingProblem.GammaScale", 10.0);
  pt.put("MySamplingProblem.GammaScale", 1.0);
  pt.put("MySamplingProblem.GammaExponent", 1.0);

  pt.put("MySamplingProblem.BetaScale", 0.0);
  pt.put("MySamplingProblem.BetaExponent", 1.0);
  
  // the target density
  auto dist = std::make_shared<Target>(1.0)->AsDensity(); 
  
  // create a sampling problem
  auto problem = std::make_shared<ExpensiveSamplingProblem>(dist, pt.get_child("MySamplingProblem"));
  
  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Zero(1);
  
  // create mcmc
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem);

  // run mcmc
  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);
  
  if( !filename.empty() ) { samps->WriteToFile(filename, dataset); }

  std::cout << samps->Covariance() << std::endl;
  std::cout << std::endl;

  return samps->RunningCovariance();
}

int main(int argc, char **argv) {
  const int ierr = MPI_Init(nullptr, nullptr);

  const unsigned int Ntrials = 10;
  const unsigned int N = 1.0e5;
  const std::string filename = "la-mcmc.h5";

  typedef std::vector<Eigen::MatrixXd> RunningCovariance;
  std::vector<RunningCovariance> runCov(Ntrials);
  runCov[0] = RunLAMCMC(N, "/trial 0", filename);
  Eigen::MatrixXd cov = runCov[0][runCov[0].size()-1];
  for( unsigned int i=1; i<runCov.size(); ++i ) {
    runCov[i] = RunLAMCMC(N, "/trial "+std::to_string(i), filename);
    assert(runCov[i].size()==runCov[0].size());
    cov = ((double)i*cov+runCov[i][runCov[i].size()-1])/(double)(i+1);
  }

  Eigen::MatrixXd errors(Ntrials, runCov[0].size());
  for( unsigned int i=0; i<runCov.size(); ++i ) {
    for( unsigned int j=0; j<runCov[i].size(); ++j ) {
      errors(i,j) = (runCov[i][j]-cov).norm();
    }
  }
  Eigen::RowVectorXd expectedError = errors.colwise().sum()/(double)Ntrials;

  auto hdf5file = std::make_shared<HDF5File>(filename);
  hdf5file->WriteMatrix("/covariance error", errors);
  hdf5file->WriteMatrix("/expected covariance error", expectedError);
  hdf5file->Close();

  MPI_Finalize();
}
