#include "MUQ/SamplingAlgorithms/MIMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    MIMCMC::MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
    : componentFactory(componentFactory),
     samples(pt.get("NumSamples",1000))
    {

      std::shared_ptr<MultiIndexSet> gridIndices = MultiIndexFactory::CreateFullTensor(componentFactory->finestIndex()->GetVector());
      const auto rootIndex = std::make_shared<MultiIndex>(gridIndices->GetMultiLength());

      // TODO: Avoid samplingProblem setup for retreiving QOI dimensions?
      Eigen::VectorXd MImean(componentFactory->samplingProblem(componentFactory->finestIndex())->blockSizesQOI.sum());
      MImean.setZero();

      for (int i = 0; i < gridIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxHighestIndex = (*gridIndices)[i];


        auto box = std::make_shared<MIMCMCBox>(componentFactory, boxHighestIndex);
        boxes.push_back(box);
        for (int samp = 0; samp < samples; samp++) {
          box->Sample();
        }
        Eigen::VectorXd sampMean = box->Mean();
        MImean += sampMean;

        std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
      }

      std::cout << "MI Mean = \n" << MImean.transpose() << std::endl;
    }

    void MIMCMC::draw() {
      std::ofstream graphfile;
      graphfile.open ("graph");
      graphfile << "digraph {" << std::endl;
      graphfile << "nodesep=1.2;" << std::endl;
      graphfile << "splines=false;" << std::endl;
      for (auto box : boxes) {
        box->draw(graphfile);
      }
      graphfile << "}" << std::endl;
      graphfile.close();
    }

  }
}
