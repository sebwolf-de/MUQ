#ifndef PARALLELMHKERNEL_H_
#define PARALLELMHKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/config.h"

// Include parcer if MUQ has MPI
#if MUQ_HAS_MPI==1
#include "parcer/Communicator.h"
#include "parcer/Queue.h"
#include "MUQ/SamplingAlgorithms/QueueWrapper.h"
#endif

namespace muq {
  namespace SamplingAlgorithms {

#if MUQ_HAS_MPI==1
    typedef parcer::Queue<std::vector<Eigen::VectorXd>, std::pair<double, Eigen::VectorXd>, QueueWrapper> QueueType;
#endif

    /// The General Parallel MCMC framework of [Calderhead 2014]
    /**
       "A general construction for parallelizing Metropolis-Hastings algorithms"

       Note, if this is run with MPI and more than one process, then all of the processes
       except rank 0 will be capture in the constructor and used for parallel density
       evaluations.
     */
    class ParallelMHKernel : public TransitionKernel {
    public:

      ParallelMHKernel(boost::property_tree::ptree const& pt,
                       std::shared_ptr<AbstractSamplingProblem> problem);

      ParallelMHKernel(boost::property_tree::ptree const& pt,
                       std::shared_ptr<AbstractSamplingProblem> problem,
                       std::shared_ptr<MCMCProposal> proposalIn);

      virtual ~ParallelMHKernel() = default;

      virtual std::shared_ptr<MCMCProposal> Proposal(){return proposal;};

      /** Generates proposals and evaluates the target density.
      */
      virtual void PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) override;

      /** Uses the results of PreStep to produce several steps of the Markov chain. */
      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual void PrintStatus(std::string prefix) const override;

      virtual double AcceptanceRate() const{return double(numAccepts)/double(numSteps);};

    protected:

      /** Computes the transition probabilities between the entries in "proposedStates".
          Assumes the target density value has already been computed and saved
          as "LogDensity" in the sample state's metadata.
      */
      Eigen::MatrixXd ComputeTransitions() const;

      /** Given a matrix of transition probabilities, this function computes the
          stationary distribution of a discrete Markov process.
      */
      Eigen::VectorXd ComputeStationary(Eigen::MatrixXd const& transProbs) const;



      std::shared_ptr<MCMCProposal> proposal;
      std::shared_ptr<MCMCProposal> zProposal;
      unsigned int propsPerStep;

      std::shared_ptr<SamplingState> zState;
      std::vector<std::shared_ptr<SamplingState>> possibleStates;

      unsigned int numSteps = 0;
      unsigned int numAccepts = 0;

#if MUQ_HAS_MPI==1
      void SetupQueue(std::shared_ptr<parcer::Communicator> comm);

      bool useParQueue;
      std::shared_ptr<QueueType> parQueue;
#endif

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
