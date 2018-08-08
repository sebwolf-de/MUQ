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

      void draw();

    private:
      std::shared_ptr<MIComponentFactory> componentFactory;
      const int numInitialSamples;
      std::vector<std::shared_ptr<MIMCMCBox>> boxes;
    };

  }
}

#endif
