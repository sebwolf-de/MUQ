#include "MUQ/SamplingAlgorithms/GMHKernel.h"

#include <parcer/Queue.h>

#include <Eigen/Eigenvalues>

#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(GMHKernel)

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : MHKernel(pt, problem),
  N(pt.get<unsigned int>("NumProposals")),
  Np1(N+1),
  M(pt.get<unsigned int>("NumAccepted", N)) {}

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) :
  MHKernel(pt, problem, proposalIn),
  N(pt.get<unsigned int>("NumProposals")),
  Np1(N+1),
  M(pt.get<unsigned int>("NumAccepted", N)) {}

GMHKernel::~GMHKernel() {}

void GMHKernel::PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) {
  // propose N steps
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;
  for( auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) {
    *it = proposal->Sample(state); 
  }

  // compute log-target
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  for( unsigned int i=0; i<Np1; ++i ) {
    R(i) = problem->LogDensity(proposedStates[i]);
    for( auto k : proposedStates ) {
      if( k==proposedStates[i] ) { continue; }
      R(i) += proposal->LogDensity(proposedStates[i], k);
    }
  }

  // compute stationary acceptance transition probability
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(Np1,Np1);
  for( unsigned int i=0; i<Np1; ++i ) {
    for( unsigned int j=0; j<Np1; ++j ) {
      if( j==i ) { continue; }
      A(i,j) = std::fmin(1.0, std::exp(R(j)-R(i)))/(double)(Np1);
      A(i,i) -= A(i,j);
    }
  }

  // compute the dominante eigen vector and then make the sum equal to 1
  const double dominateEigenval = PowerIteration(A.transpose());
  assert(std::fabs(dominateEigenval-1.0)<1.0e-10);
  stationaryAcceptance /= stationaryAcceptance.sum();

  // compute the cumulative sum
  for( unsigned int i=1; i<Np1; ++i ) { stationaryAcceptance(i) += stationaryAcceptance(i-1); }
}

double GMHKernel::PowerIteration(Eigen::MatrixXd const& A) {
  stationaryAcceptance = Eigen::VectorXd::Ones(A.cols()).normalized();
  
  int counter = 0;
  double error=100.0*tol;
  while( error>tol && counter++<maxIt) {
    Eigen::VectorXd temp = stationaryAcceptance;
    stationaryAcceptance = (A*temp).normalized();
    error = (temp-stationaryAcceptance).stableNorm();
  }
  
  // return the dominate eigenvalue
  return stationaryAcceptance.transpose()*A*stationaryAcceptance;
}

std::vector<std::shared_ptr<SamplingState> > GMHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> state) {
  std::vector<std::shared_ptr<SamplingState> > newStates(M, nullptr);

  // sample the new states
  for( auto it=newStates.begin(); it!=newStates.end(); ++it ) {
    // determine which proposed state to return
    const double uniform = RandomGenerator::GetUniform();
    unsigned int index = 0;
    for( ; index<Np1; ++index ) { if( uniform<stationaryAcceptance(index) ) { break; }  }

    // store it
    *it = proposedStates[index];
  }

  return newStates;
}

Eigen::VectorXd GMHKernel::CumulativeStationaryAcceptance() const {
  // make sure this object has been populated
  assert(stationaryAcceptance.size()==Np1);
  
  return stationaryAcceptance;
}
