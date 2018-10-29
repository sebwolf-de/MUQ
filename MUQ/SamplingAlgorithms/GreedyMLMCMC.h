#ifndef GreedyMLMCMC_H
#define GreedyMLMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    class GreedyMLMCMC {
    public:
      GreedyMLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      void run();

      Eigen::VectorXd meanQOI();

      void draw(bool drawSamples = true);

      std::shared_ptr<MIMCMCBox> getBox(int index);

    private:
      const double e;
      const double beta;
      const int levels;
      std::shared_ptr<MIComponentFactory> componentFactory;
      const int numInitialSamples;
      std::vector<std::shared_ptr<MIMCMCBox>> boxes;
    };

  }
}

#endif
