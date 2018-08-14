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
  M(pt.get<unsigned int>("NumAccpeted", N)) {}

GMHKernel::GMHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) :
  MHKernel(pt, problem, proposalIn),
  N(pt.get<unsigned int>("NumProposals")),
  Np1(N+1),
  M(pt.get<unsigned int>("NumAccpeted", N)) {}

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
      A(i,j) = std::fmin(1.0, R(j)/R(i))/(double)(Np1);
      A(i,i) -= A(i,j);
    }
  }

  Eigen::EigenSolver<Eigen::MatrixXd> eig(A);
  // the larget eigenvalue should be 1
  assert(std::fabs(eig.eigenvalues() (0).real()-1.0)<1.0e-10);
  assert(std::fabs(eig.eigenvalues() (0).imag())<1.0e-10);

  stationaryAcceptance = eig.eigenvectors().col(0).real();
  stationaryAcceptance /= stationaryAcceptance.sum();
  
  // compute the cumulative sum
  for( unsigned int i=1; i<Np1; ++i ) { stationaryAcceptance(i) += stationaryAcceptance(i-1); }
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
