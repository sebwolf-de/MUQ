#ifndef MIMCMC_H
#define MIMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Multiindex MCMC method.
        @details A basic MIMCMC method based on a fixed
        number of samples for all model indices.
     */
    class MIMCMC : public SamplingAlgorithm {
    public:
      MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      virtual std::shared_ptr<SampleCollection> RunImpl() override;

      virtual std::shared_ptr<SampleCollection> GetSamples() const;
      virtual std::shared_ptr<SampleCollection> GetQOIs() const;

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
