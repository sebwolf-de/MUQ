#include "MUQ/OptimalExperimentalDesign/LogDifference.h"

using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

LogDifference::LogDifference(std::shared_ptr<Distribution> const& like, std::shared_ptr<Distribution> const& evidence) :  ModPiece(Eigen::VectorXi::Ones(3), Eigen::VectorXi::Ones(1)), likelihood(like), evidence(evidence) {}

void LogDifference::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  const Eigen::VectorXd& x = inputs[0];
  const Eigen::VectorXd& y = inputs[1];
  const Eigen::VectorXd& d = inputs[2];

  outputs.resize(1);
  outputs[0] = Eigen::VectorXd::Constant(1, likelihood->LogDensity(y, x, d)-evidence->LogDensity(y, d));
}
