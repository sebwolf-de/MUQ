#include "MUQ/OptimalExperimentalDesign/LogDifference.h"

using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

LogDifference::LogDifference(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence) :  ModPiece(Eigen::Vector3i(evidence->varSize,likelihood->hyperSizes(0), evidence->hyperSizes(0)), Eigen::VectorXi::Ones(1)), likelihood(likelihood), evidence(evidence) {}

void LogDifference::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  const Eigen::VectorXd& y = inputs[0];
  const Eigen::VectorXd& x = inputs[1];
  const Eigen::VectorXd& d = inputs[2];

  outputs.resize(1);
  outputs[0] = Eigen::VectorXd::Constant(1, likelihood->LogDensity(y, x, d)-evidence->LogDensity(y, d));
}
