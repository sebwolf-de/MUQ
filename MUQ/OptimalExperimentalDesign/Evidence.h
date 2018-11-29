#ifndef EVIDENCE_H_
#define EVIDENCE_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/config.h"

#if MUQ_HAS_PARCER==1
#include <parcer/Communicator.h>
#endif

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace OptimalExperimentalDesign {

    class Evidence : public muq::Modeling::Distribution {
    public:

      Evidence(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& biasing, boost::property_tree::ptree pt);

#if MUQ_HAS_PARCER==1
      Evidence(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& biasing, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> const& comm);
#endif

      virtual ~Evidence() = default;
    private:

      void CreateGraph(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood);

      virtual double LogDensityImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      const unsigned int numImportanceSamples;

      std::shared_ptr<muq::Modeling::Distribution> biasing;

      std::shared_ptr<muq::Modeling::WorkGraph> graph;

#if MUQ_HAS_PARCER==1
      std::shared_ptr<parcer::Communicator> comm = nullptr;
#endif

    };
  } // namespace OptimalExperimentalDesign
} // namespace muq

#endif
