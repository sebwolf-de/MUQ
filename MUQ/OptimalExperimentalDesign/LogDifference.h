#ifndef LOGDIFFERENCE_H_
#define LOGDIFFERENCE_H_

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace OptimalExperimentalDesign {
    class LogDifference : public muq::Modeling::ModPiece {
    public:
      LogDifference(std::shared_ptr<muq::Modeling::Distribution> const& evidence);


      virtual ~LogDifference() = default;
    private:
      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      std::shared_ptr<muq::Modeling::Distribution> evidence;
    };
  } // namespace OptimalExperimentalDesign
} // namespace muq

#endif
