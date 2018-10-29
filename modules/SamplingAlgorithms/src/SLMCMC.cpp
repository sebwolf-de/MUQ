#include "MUQ/SamplingAlgorithms/SLMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    SLMCMC::SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
     : componentFactory(componentFactory)
    {
      auto index = componentFactory->finestIndex();

      auto problem = componentFactory->samplingProblem(index);
      auto proposal = componentFactory->proposal(index, problem);

      pt::ptree ptBlockID;
      ptBlockID.put("BlockIndex",0);
      std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
      kernels[0] = std::make_shared<MHKernel>(ptBlockID,problem,proposal);

      Eigen::VectorXd startingPoint = componentFactory->startingPoint(index);

      coarse_chain = std::make_shared<SingleChainMCMC>(pt,kernels,startingPoint);
    }

    void SLMCMC::run() {
      coarse_chain->Run();
    }

    Eigen::VectorXd SLMCMC::meanQOI() {
      return coarse_chain->GetQOIs()->Mean();
    }

    Eigen::VectorXd SLMCMC::meanParameter() {
      auto samps = coarse_chain->GetSamples();
      return samps->Mean();
    }

  }
}
