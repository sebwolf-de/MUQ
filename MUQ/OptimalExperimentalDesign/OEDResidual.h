#ifndef OEDRESIDUAL_H_
#define OEDRESIDUAL_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/config.h"

#if MUQ_HAS_PARCER==1
#include <parcer/Communicator.h>
#endif

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/Distributions/Distribution.h"

#include "MUQ/Approximation/Regression/LocalRegression.h"

namespace muq {
  namespace OptimalExperimentalDesign {
    class OEDResidual : public muq::Modeling::ModPiece {
    public:

      /// Use the likelihood as the baising distribution---Monte Carlo estimate
      OEDResidual(std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& evidence, boost::property_tree::ptree pt);

#if MUQ_HAS_PARCER==1
      /// Use the likelihood as the baising distribution---Monte Carlo estimate
      OEDResidual(std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& evidence, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> const& comm);
#endif

      ~OEDResidual() = default;
    private:

      void CreateGraph();

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      void RandomlyRefineNear(std::shared_ptr<muq::Approximation::LocalRegression> const& reg, Eigen::VectorXd const& xd, double const radius);

      void RefineAt(std::shared_ptr<muq::Approximation::LocalRegression> const& reg, Eigen::VectorXd const& pnt, double const radius);

      const unsigned int numImportanceSamples;

      std::shared_ptr<muq::Modeling::Distribution> likelihood;
      std::shared_ptr<muq::Modeling::Distribution> evidence;

      unsigned int totalRefinements = 0;

      const bool bruteForce;

      const double gamma0;
      const double radius0;

      std::shared_ptr<muq::Modeling::WorkGraph> graph;

      boost::property_tree::ptree pt_regression;

#if MUQ_HAS_PARCER==1
      std::shared_ptr<parcer::Communicator> comm = nullptr;
#endif
    };
  } // namespace OptimalExperimentalDesign
} // namespace muq

#endif
