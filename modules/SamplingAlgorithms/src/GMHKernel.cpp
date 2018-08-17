#include "MUQ/SamplingAlgorithms/GMHKernel.h"

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

void GMHKernel::SerialProposal(std::shared_ptr<SamplingState> state) {
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;
  for( auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) { *it = proposal->Sample(state); }
}

#if MUQ_HAS_PARCER
void GMHKernel::ParallelProposal(std::shared_ptr<SamplingState> state) {
  // if we only have one processor, just propose in serial
  if( comm->GetSize()==1 ) { return SerialProposal(state); }

  auto proposalQueue = std::make_shared<ProposalQueue>(proposal, comm);

  if( comm->GetRank()==0 ) {
    assert(state);
    
    std::vector<unsigned int> proposalIDs(N);
    for( auto id=proposalIDs.begin(); id!=proposalIDs.end(); ++id ) { *id = proposalQueue->SubmitWork(state); }

    proposedStates.resize(Np1, nullptr);
    proposedStates[0] = state;
    for( unsigned int i=0; i<N; ++i ) {
      proposedStates[i+1] = proposalQueue->GetResult(proposalIDs[i]);
    }
  }
}
#endif

void GMHKernel::SerialAcceptanceDensity() {
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  for( unsigned int i=0; i<Np1; ++i ) {
    R(i) = problem->LogDensity(proposedStates[i]);
    for( auto k : proposedStates ) {
      if( k==proposedStates[i] ) { continue; }
      R(i) += proposal->LogDensity(proposedStates[i], k);
    }
  }

  CumulativeAcceptanceDensity(R);
}

#if MUQ_HAS_PARCER
void GMHKernel::ParallelLogTarget(Eigen::VectorXd& R) {
  auto logTargetQueue = std::make_shared<LogTargetQueue>(problem, comm);

  if( comm->GetRank()==0 ) {
    std::vector<unsigned int> logTargetIDs(Np1);
    for( unsigned int i=0; i<Np1; ++i ) { logTargetIDs[i] = logTargetQueue->SubmitWork(proposedStates[i]); }

    R = Eigen::VectorXd::Zero(Np1);
    for( unsigned int i=0; i<Np1; ++i ) { R(i) = logTargetQueue->GetResult(logTargetIDs[i]);}
  }
}

struct EvaluateProposalDensity {
  inline EvaluateProposalDensity(std::shared_ptr<MCMCProposal> proposal) : proposal(proposal) {}

  inline double Evaluate(std::pair<std::shared_ptr<SamplingState>, std::shared_ptr<SamplingState> > states) {
    return proposal->LogDensity(states.first, states.second);
  }

  std::shared_ptr<MCMCProposal> proposal;
};

void GMHKernel::ParallelAcceptanceDensity() {
  // if we only have one processor, compute log-target in serial
  if( comm->GetSize()==1 ) { return SerialAcceptanceDensity(); }

  /*Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  if( comm->GetRank()==0 ) {
    for( unsigned int i=0; i<Np1; ++i ) { R(i) = problem->LogDensity(proposedStates[i]); }
  }*/
  Eigen::VectorXd R;
  ParallelLogTarget(R);

  auto evalPropDens = std::make_shared<EvaluateProposalDensity>(proposal);
  parcer::Queue<std::pair<std::shared_ptr<SamplingState>, std::shared_ptr<SamplingState> >, double, EvaluateProposalDensity> proposalDensityQueue(evalPropDens, comm);

  if( comm->GetRank()==0 ) {
    std::vector<unsigned int> ids(Np1*N);
    unsigned int cnt = 0;
    for( unsigned int i=0; i<Np1; ++i ) {
      for( auto k : proposedStates ) {
	if( k==proposedStates[i] ) { continue; }
	ids[cnt++] = proposalDensityQueue.SubmitWork(std::pair<std::shared_ptr<SamplingState>, std::shared_ptr<SamplingState> >(proposedStates[i], k));
      }
    }

    cnt = 0;
    for( unsigned int i=0; i<Np1; ++i ) {
      for( unsigned int k=0; k<N; ++k ) {
	if( k==i ) { continue; }
	R(i) += proposalDensityQueue.GetResult(ids[cnt++]);
      }
    }

    /*for( unsigned int i=0; i<Np1; ++i ) {
      for( auto k : proposedStates ) {
	if( k==proposedStates[i] ) { continue; }
	R(i) += proposal->LogDensity(proposedStates[i], k);
      }
      }*/
    
    CumulativeAcceptanceDensity(R);
  }
}
#endif

Eigen::MatrixXd GMHKernel::AcceptanceMatrix(Eigen::VectorXd const& R) const {
  // compute stationary acceptance transition probability
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(Np1,Np1);
  for( unsigned int i=0; i<Np1; ++i ) {
    for( unsigned int j=0; j<Np1; ++j ) {
      if( j==i ) { continue; }
      A(i,j) = std::fmin(1.0, std::exp(R(j)-R(i)))/(double)(Np1);
      A(i,i) -= A(i,j);
    }
  }

  return A;
}

void GMHKernel::CumulativeAcceptanceDensity(Eigen::VectorXd const& R) {
  const Eigen::MatrixXd& A = AcceptanceMatrix(R);
  
  // compute the dominante eigen vector and then make the sum equal to 1
  const double dominateEigenval = PowerIteration(A.transpose());
  assert(std::fabs(dominateEigenval-1.0)<1.0e-10);
  stationaryAcceptance /= stationaryAcceptance.sum();
  
  // compute the cumulative sum
  for( unsigned int i=1; i<Np1; ++i ) { stationaryAcceptance(i) += stationaryAcceptance(i-1); }
}

void GMHKernel::PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) {
  // propose N steps
  if( comm ) {
#if MUQ_HAS_PARCER
    if( comm->GetRank()==0 ) { assert(state); }
    ParallelProposal(state);
#else
    SerialProposal(state);
#endif
  } else { SerialProposal(state); }

  // compute cumulative acceptance density
  if( comm ) { 
#if MUQ_HAS_PARCER
    ParallelAcceptanceDensity();
#else
    SerialAcceptanceDensity();
#endif
  } else { SerialAcceptanceDensity(); }
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

std::vector<std::shared_ptr<SamplingState> > GMHKernel::SampleStationary(std::shared_ptr<SamplingState> state) {
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

std::vector<std::shared_ptr<SamplingState> > GMHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> state) {
  if( !comm ) { return SampleStationary(state); } 
  
#if MUQ_HAS_PARCER
  return comm->GetRank()==0 ? SampleStationary(state) : std::vector<std::shared_ptr<SamplingState> >(M, nullptr);
#else
  return SampleStationary(state);
#endif
}

Eigen::VectorXd GMHKernel::CumulativeStationaryAcceptance() const {
  // make sure this object has been populated
  assert(stationaryAcceptance.size()==Np1);
  
  return stationaryAcceptance;
}

#if MUQ_HAS_PARCER
std::shared_ptr<parcer::Communicator> GMHKernel::GetCommunicator() const { return comm; }
#endif
