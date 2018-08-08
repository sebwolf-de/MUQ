#ifndef SLMCMC_H
#define SLMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    class SLMCMC {

    public:
      SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

    private:
      std::shared_ptr<MIComponentFactory> componentFactory;

    };

  }
}

#endif
