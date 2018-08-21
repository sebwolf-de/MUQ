#include "MUQ/SamplingAlgorithms/GMHKernel.h"

#include <Eigen/Eigenvalues>

#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(GMHKernel)

#if MUQ_HAS_PARCER
typedef std::pair<std::shared_ptr<SamplingState>, bool> CurrentState;

struct ProposeState {
  inline ProposeState(std::shared_ptr<MCMCProposal> proposal, std::shared_ptr<AbstractSamplingProblem> problem) : proposal(proposal), problem(problem) {}

  inline std::shared_ptr<SamplingState> Evaluate(CurrentState state) {
    std::shared_ptr<SamplingState> proposed = state.second ? proposal->Sample(state.first) : state.first;
    proposed->meta["log-target"] = problem->LogDensity(proposed);
    return proposed;
  }

  std::shared_ptr<MCMCProposal> proposal;
  std::shared_ptr<AbstractSamplingProblem> problem;
};

typedef parcer::Queue<CurrentState, std::shared_ptr<SamplingState>, ProposeState> ProposalQueue;
#endif

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
  // propose the points
  proposedStates.resize(Np1, nullptr);
  proposedStates[0] = state;
  proposedStates[0]->meta["log-target"] = problem->LogDensity(state);
  for( auto it = proposedStates.begin()+1; it!=proposedStates.end(); ++it ) {
    *it = proposal->Sample(state);
    (*it)->meta["log-target"] = problem->LogDensity(*it);
  }

  // evaluate the target density
  Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
  for( unsigned int i=0; i<Np1; ++i ) { R(i) = boost::any_cast<double const>(proposedStates[i]->meta["log-target"]); }

  // compute stationary transition probability
  AcceptanceDensity(R);
}

#if MUQ_HAS_PARCER
void GMHKernel::ParallelProposal(std::shared_ptr<SamplingState> state) {
  // if we only have one processor, just propose in serial
  if( comm->GetSize()==1 ) { return SerialProposal(state); }

  // create a queue to propose and evaluate the log-target
  auto helper = std::make_shared<ProposeState>(proposal, problem);
  auto proposalQueue = std::make_shared<ProposalQueue>(helper, comm);
  
  if( comm->GetRank()==0 ) {
    assert(state);

    // submit the work
    std::vector<unsigned int> proposalIDs(Np1);
    proposalIDs[0] = proposalQueue->SubmitWork(CurrentState(state, false)); // evaluate the current state
    for( auto id=proposalIDs.begin()+1; id!=proposalIDs.end(); ++id ) { *id = proposalQueue->SubmitWork(CurrentState(state, true)); }

    // retrieve the work
    proposedStates.resize(Np1, nullptr);
    Eigen::VectorXd R = Eigen::VectorXd::Zero(Np1);
    for( unsigned int i=0; i<Np1; ++i ) {
      std::shared_ptr<SamplingState> evalState = proposalQueue->GetResult(proposalIDs[i]);
      
      proposedStates[i] = evalState;
      R(i) = boost::any_cast<double const>(evalState->meta["log-target"]);
    }

    // compute stationary transition probability
    AcceptanceDensity(R);
  }
}
#endif
  
void GMHKernel::AcceptanceDensity(Eigen::VectorXd& R) {
  // update log-target with proposal density
  for( unsigned int i=0; i<Np1; ++i ) {
    for( auto k : proposedStates ) {
      if( k==proposedStates[i] ) { continue; }
      R(i) += proposal->LogDensity(proposedStates[i], k);
    }
  }

  // compute the cumlative acceptance density
  ComputeStationaryAcceptance(R);
}

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

void GMHKernel::ComputeStationaryAcceptance(Eigen::VectorXd const& R) {
  const Eigen::MatrixXd& A = AcceptanceMatrix(R);

  stationaryAcceptance = Eigen::VectorXd::Ones(A.cols()).normalized();

  Eigen::MatrixXd mat(Np1+1, Np1);
  mat.block(0,0,Np1,Np1) = A.transpose()-Eigen::MatrixXd::Identity(Np1,Np1);
  mat.row(Np1) = Eigen::RowVectorXd::Ones(Np1);

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(Np1+1);
  rhs(Np1) = 1.0;

  stationaryAcceptance = mat.colPivHouseholderQr().solve(rhs);
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
}

std::vector<std::shared_ptr<SamplingState> > GMHKernel::SampleStationary() const {
  std::vector<std::shared_ptr<SamplingState> > newStates(M, nullptr);

  // get the indices of the proposed states that are accepted
  Eigen::VectorXi indices = RandomGenerator::GetDiscrete(stationaryAcceptance, M);
  assert(indices.size()==M);

  // store the accepted states
  for( unsigned int i=0; i<M; ++i ) { newStates[i] = proposedStates[indices[i]]; }

  return newStates;
}

std::vector<std::shared_ptr<SamplingState> > GMHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> state) {
  if( !comm ) { return SampleStationary(); } 
  
#if MUQ_HAS_PARCER
  return comm->GetRank()==0 ? SampleStationary() : std::vector<std::shared_ptr<SamplingState> >(M, nullptr);
#else
  return SampleStationary();
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
