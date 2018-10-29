#ifndef MIMCMC_H
#define MIMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    class MIMCMC {
    public:
      MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      void run();

      Eigen::VectorXd meanQOI();

      void draw(bool drawSamples = true);

    private:
      std::shared_ptr<MultiIndexSet> gridIndices;
      std::shared_ptr<MIComponentFactory> componentFactory;
      const int samples;
      std::vector<std::shared_ptr<MIMCMCBox>> boxes;

    };

  }
}

#endif
