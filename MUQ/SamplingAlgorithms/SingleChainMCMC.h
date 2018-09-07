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

      std::shared_ptr<SamplingState> SaveSamples(std::vector<std::shared_ptr<SamplingState> > const& newStates, unsigned int& sampNum) const;

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
