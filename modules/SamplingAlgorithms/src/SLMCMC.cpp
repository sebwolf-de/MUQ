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

      auto coarse_chain = std::make_shared<SingleChainMCMC>(pt,kernels,startingPoint);
      coarse_chain->Run();
      auto samps = coarse_chain->GetQOIs();

      Eigen::VectorXd sampMean = samps->Mean();

      std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
    }

  }
}
