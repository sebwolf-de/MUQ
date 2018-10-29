#include "MUQ/SamplingAlgorithms/MIMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    MIMCMC::MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
    : componentFactory(componentFactory),
     samples(pt.get("NumSamples",1000))
    {

      gridIndices = MultiIndexFactory::CreateFullTensor(componentFactory->finestIndex()->GetVector());

      for (int i = 0; i < gridIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxHighestIndex = (*gridIndices)[i];

        auto box = std::make_shared<MIMCMCBox>(componentFactory, boxHighestIndex);
        boxes.push_back(box);
      }

    }

    void MIMCMC::run() {

      for (auto box : boxes) {
        for (int samp = 0; samp < samples; samp++) {
          box->Sample();
        }
      }

    }

    Eigen::VectorXd MIMCMC::meanQOI() {
      // TODO: Avoid samplingProblem setup for retreiving QOI dimensions?
      // Compute full QOI estimate
      Eigen::VectorXd MImean(componentFactory->samplingProblem(componentFactory->finestIndex())->blockSizesQOI.sum());
      MImean.setZero();

      for (auto box : boxes) {
        Eigen::VectorXd sampMean = box->MeanQOI();

        MImean += sampMean;
      }

      return MImean;
    }

    void MIMCMC::draw(bool drawSamples) {
      std::ofstream graphfile;
      graphfile.open ("graph");
      graphfile << "digraph {" << std::endl;
      graphfile << "nodesep=1.2;" << std::endl;
      graphfile << "splines=false;" << std::endl;
      for (auto box : boxes) {
        box->draw(graphfile, drawSamples);
      }
      graphfile << "}" << std::endl;
      graphfile.close();
    }

  }
}
