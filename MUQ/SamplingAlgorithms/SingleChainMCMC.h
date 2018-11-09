#ifndef SINGLECHAINMCMC_H
#define SINGLECHAINMCMC_H

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"
#include "MUQ/SamplingAlgorithms/TransitionKernel.h"
#include "MUQ/SamplingAlgorithms/ThinScheduler.h"

#include <vector>

#include <boost/property_tree/ptree.hpp>

namespace muq{
  namespace SamplingAlgorithms{

    /** @ingroup MCMC
        @class SingleChainMCMC
        @brief Defines an MCMC sampler with a single chain
        @details 
        <B>Configuration Parameters:</B>
        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "NumSamples"  | Int | - | The total number of steps (including burnin) to take, i.e., the length of the Markov chain. |
        "BurnIn"      | Int | 0 | The number of steps at the beginning of the chain to ignore. |
        "PrintLevel"  | Int | 3 | The amount of information to print to std::cout. Valid values are in [0,1,2,3] with  0 = Nothing, 3 = The most |
        "KernelList"  | String | - | A comma separated list of other parameter blocks that define the transition kernels for each Metropolis-in-Gibbs block. |
    */
    class SingleChainMCMC : public SamplingAlgorithm
    {

    public:

      SingleChainMCMC(boost::property_tree::ptree pt, std::shared_ptr<AbstractSamplingProblem> problem);

#if MUQ_HAS_PARCER
      SingleChainMCMC(boost::property_tree::ptree pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<parcer::Communicator> comm);
#endif

      SingleChainMCMC(boost::property_tree::ptree              pt,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      std::vector<std::shared_ptr<TransitionKernel>> kernelsIn);

      virtual ~SingleChainMCMC() = default;

      virtual std::vector<std::shared_ptr<TransitionKernel>>& Kernels(){return kernels;};

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

    protected:

      std::shared_ptr<SamplingState> SaveSamples(std::vector<std::shared_ptr<SamplingState> > const& newStates, std::shared_ptr<SamplingState>& lastSavedState, unsigned int& sampNum) const;

      bool ShouldSave(unsigned int const sampNum) const;

      void PrintStatus(unsigned int currInd) const{PrintStatus("",currInd);};
      void PrintStatus(std::string prefix, unsigned int currInd) const;

      std::shared_ptr<SaveSchedulerBase> scheduler;

      unsigned int numSamps;
      unsigned int burnIn;
      unsigned int printLevel;

      // A vector of transition kernels: One for each block
      std::vector<std::shared_ptr<TransitionKernel>> kernels;

    private:

      void SetUp(boost::property_tree::ptree pt, std::shared_ptr<AbstractSamplingProblem> problem);
    }; // class SingleChainMCMC
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // #ifndef SINGLECHAINMCMC_H
