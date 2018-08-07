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

      SingleChainMCMC(boost::property_tree::ptree              pt,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      std::vector<Eigen::VectorXd> x0);

      SingleChainMCMC(boost::property_tree::ptree              pt,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      Eigen::VectorXd x0)
                    : SingleChainMCMC(pt, problem, std::vector<Eigen::VectorXd>(1,x0)) {};

      SingleChainMCMC(boost::property_tree::ptree& pt,
                      std::vector<std::shared_ptr<TransitionKernel>> kernels,
                      std::vector<Eigen::VectorXd> x0);

      SingleChainMCMC(boost::property_tree::ptree& pt,
                      std::vector<std::shared_ptr<TransitionKernel>> kernels,
                      Eigen::VectorXd x0)
                    : SingleChainMCMC(pt, kernels, std::vector<Eigen::VectorXd>(1,x0)) {};


      virtual ~SingleChainMCMC() = default;

      virtual std::vector<std::shared_ptr<TransitionKernel>>& Kernels(){return kernels;};

      virtual std::shared_ptr<SampleCollection> RunImpl() override;

      virtual void Sample();

      virtual double TotalTime() { return totalTime; }

    protected:

      std::shared_ptr<SaveSchedulerBase> scheduler;
      std::shared_ptr<SaveSchedulerBase> schedulerQOI;

      void PrintStatus(unsigned int currInd) const{PrintStatus("",currInd);};
      void PrintStatus(std::string prefix, unsigned int currInd) const;

      unsigned int numSamps;
      unsigned int burnIn;
      unsigned int printLevel;

      // A vector of transition kernels: One for each block
      std::vector<std::shared_ptr<TransitionKernel>> kernels;

    private:

      unsigned int sampNum = 1;
      std::shared_ptr<SamplingState> prevState = nullptr;
      std::shared_ptr<SamplingState> lastSavedState = nullptr;
      std::shared_ptr<SamplingState> lastSavedQOI = nullptr;
      double totalTime = 0.0;

    }; // class SingleChainMCMC

  } // namespace SamplingAlgorithms
} // namespace muq

#endif // #ifndef SINGLECHAINMCMC_H
