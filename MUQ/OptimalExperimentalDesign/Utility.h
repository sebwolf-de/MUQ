#ifndef UTILITY_H_
#define UTILITY_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/config.h"

#if MUQ_HAS_PARCER==1
#include <parcer/Communicator.h>
#endif

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace OptimalExperimentalDesign {
    class Utility : public muq::Modeling::ModPiece {
    public:
      Utility(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& evidence, std::shared_ptr<muq::Modeling::Distribution> const& biasing, boost::property_tree::ptree pt);

#if MUQ_HAS_PARCER==1
      Utility(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& evidence, std::shared_ptr<muq::Modeling::Distribution> const& biasing, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> const& comm);
#endif

      virtual ~Utility() = default;
    private:
      void CreateGraph(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& evidence);

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      const unsigned int numImportanceSamples;

      std::shared_ptr<muq::Modeling::Distribution> evidence;

      std::shared_ptr<muq::Modeling::Distribution> biasing;

#if MUQ_HAS_PARCER==1
      std::shared_ptr<parcer::Communicator> comm;
#endif

      std::shared_ptr<muq::Modeling::WorkGraph> graph;
    };
  } // namespace OptimalExperimentalDesign
} // namespace muq

#endif
