#include "MUQ/OptimalExperimentalDesign/LogDifference.h"

using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

LogDifference::LogDifference(std::shared_ptr<Distribution> const& evidence) :  ModPiece(Eigen::Vector3i(evidence->varSize, evidence->hyperSizes(0), 1), Eigen::VectorXi::Ones(1)), evidence(evidence) {}

void LogDifference::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  const Eigen::VectorXd& y = inputs[0];
  const Eigen::VectorXd& d = inputs[1];
  const double loglike = inputs[2] (0);

  outputs.resize(1);
  outputs[0] = Eigen::VectorXd::Constant(1, loglike-evidence->LogDensity(y, d));
}
