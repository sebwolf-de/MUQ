#include "MUQ/SamplingAlgorithms/ParallelMHKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <Eigen/Dense>

#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(ParallelMHKernel)

ParallelMHKernel::ParallelMHKernel(pt::ptree const& pt,
                                   std::shared_ptr<AbstractSamplingProblem> problem) :
                                   TransitionKernel(pt, problem) {

  // Extract the proposal parts from the ptree
  std::string proposalName = pt.get<std::string>("Proposal");
  std::string zName = pt.get<std::string>("Auxillary Proposal");

  boost::property_tree::ptree subTree = pt.get_child(proposalName);
  subTree.put("BlockIndex", blockInd);

  // Construct the proposal
  proposal = MCMCProposal::Construct(subTree, problem);

  subTree = pt.get_child(zName);
  subTree.put("BlockIndex", blockInd);
  zProposal = MCMCProposal::Construct(subTree, problem);

  propsPerStep = pt.get("ProposalsPerStep", 1);

#if MUQ_HAS_MPI==1
  // If we want to evaluate the posterior with MPI parallelism, we need to set up the queue
  if(pt.get("ParallelEvaluations",false)){
    auto comm = std::make_shared<parcer::Communicator>();
    SetupQueue(comm);
  }
#endif

}

ParallelMHKernel::ParallelMHKernel(pt::ptree const& pt,
                                   std::shared_ptr<AbstractSamplingProblem> problem,
                                   std::shared_ptr<MCMCProposal> proposalIn) :
                                   TransitionKernel(pt, problem), proposal(proposalIn)
{
  propsPerStep = pt.get("ProposalsPerStep", 1);

  // If we want to evaluate the posterior with MPI parallelism, we need to set up the queue
  if(pt.get("ParallelEvaluations",false)){
    SetupQueue();
  }
}



void ParallelMHKernel::PreStep(unsigned int const t, std::shared_ptr<SamplingState> prevState)
{
  possibleStates.resize(propsPerStep+1);

  // zState = zProposal->Sample(prevState);
  // zState->meta["LogTarget"] = problem->LogDensity(prevState);

  // The first possible state is the current state
  if(!prevState->HasMeta("LogTarget") )
    prevState->meta["LogTarget"] = problem->LogDensity(prevState);
  possibleStates.at(0) = prevState;

  // Propose a bunch of points
  for(int i=1; i<propsPerStep+1; ++i)
    possibleStates.at(i) = proposal->Sample(prevState);

  // Compute the target density for each state
  for(int i=1; i<propsPerStep+1; ++i)
    possibleStates.at(i)->meta["LogTarget"] = problem->LogDensity(possibleStates.at(i));
}

void ParallelMHKernel::PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state){
  proposal->Adapt(t,state);
}

std::vector<std::shared_ptr<SamplingState>> ParallelMHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> prevState){

  assert(proposal);

  // Compute the stationary probability distribution over the proposed states
  ComputeTransitions();

  Eigen::VectorXd stationaryProbs = ComputeStationary(ComputeTransitions());
  //std::cout << "Stationary probabilities = \n" << stationaryProbs.transpose() << std::endl;
  Eigen::VectorXi propInds = RandomGenerator::GetDiscrete(stationaryProbs, propsPerStep);

  std::vector<std::shared_ptr<SamplingState>> outputs(propsPerStep);

  for(int i=0; i<propsPerStep; ++i){
    outputs.at(i) = possibleStates.at(propInds(i));

    numSteps++;
    if(i>0){
      numAccepts += outputs.at(i) != outputs.at(i-1);
    }else{
      numAccepts += outputs.at(i) != prevState;
    }
  }

  return outputs;
}

void ParallelMHKernel::PrintStatus(const std::string prefix) const
{
  std::stringstream msg;
  msg << std::setprecision(2);
  msg << prefix << "Acceptance Rate = "  << 100.0*double(numAccepts)/double(numSteps) << "%";

  std::cout << msg.str() << std::endl;
}

Eigen::MatrixXd ParallelMHKernel::ComputeTransitions() const
{
  const unsigned int numStates = possibleStates.size();

  Eigen::MatrixXd transProbs = Eigen::MatrixXd::Zero(numStates, numStates);

  Eigen::VectorXd allLogParts = Eigen::VectorXd::Zero(numStates);
  for(int j=0; j<numStates; ++j){
    assert(possibleStates.at(j));
    assert(possibleStates.at(j)->HasMeta("LogTarget"));

    for(int notj = 0; notj<j; ++notj)
      allLogParts(j) += proposal->LogDensity(possibleStates.at(j),possibleStates.at(notj));
    for(int notj = j+1; notj<numStates; ++notj)
      allLogParts(j) += proposal->LogDensity(possibleStates.at(j),possibleStates.at(notj));
  }

  for(int j=0; j<numStates; ++j){
    double pj = boost::any_cast<double>( possibleStates.at(j)->meta["LogTarget"] );

    for(int i=0; i<numStates; ++i){

      if(i!=j){
        double pi = boost::any_cast<double>( possibleStates.at(i)->meta["LogTarget"] );

        transProbs(i,j) = (1.0/(numStates-1.0)) * std::min(1.0, std::exp( pj - pi + allLogParts(j) - allLogParts(i)));
      }
    }
  }

  for(int i=0; i<numStates; ++i){
    double rowSum = transProbs.row(i).sum();
    if(rowSum > 1.0){
      transProbs /= rowSum;
      rowSum = 1.0;
    }
    transProbs(i,i) = 1.0 - rowSum;
    assert(transProbs(i,i)>=0);
  }

  return transProbs;
}


Eigen::VectorXd ParallelMHKernel::ComputeStationary(Eigen::MatrixXd const& transProbs) const
{
  int numPossible = transProbs.rows();

  std::cout << 1.0 - transProbs.diagonal().mean() << std::endl;

  Eigen::MatrixXd A(numPossible+1,numPossible);
  A.block(0,0,numPossible,numPossible) = transProbs.transpose() - Eigen::MatrixXd::Identity(numPossible,numPossible);
  A.row(numPossible) = Eigen::RowVectorXd::Ones(numPossible);

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(numPossible+1);
  rhs(numPossible) = 1.0;

  return A.colPivHouseholderQr().solve(rhs);

}

void ParallelMHKernel::SetupQueue(std::shared_ptr<parcer::Communicator> comm)
{
#if MUQ_HAS_MPI==1
  useParQueue = true;
  parQueue = std::make_shared<QueueType>(, comm);
#else

#endif
}
