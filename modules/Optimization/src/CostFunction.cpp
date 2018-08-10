#include "MUQ/Optimization/CostFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

CostFunction::CostFunction(Eigen::VectorXi const& inputSizes) : ModPiece(inputSizes, Eigen::VectorXi::Ones(1)) {}

CostFunction::~CostFunction() {}

void CostFunction::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  outputs.resize(1);
  outputs.at(0) = Eigen::VectorXd::Constant(1, CostImpl(input));
}

void CostFunction::GradientImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  GradientImpl(inputDimWrt, input, sensitivity);
}

void CostFunction::GradientImpl(unsigned int const inputDimWrt, ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  ModPiece::GradientImpl(0, inputDimWrt, input, sensitivity);
}

double CostFunction::Cost(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) {
  return Evaluate(input).at(0) (0);
}

Eigen::VectorXd const& CostFunction::Gradient(unsigned int const inputDimWrt, std::vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  return ModPiece::Gradient(0, inputDimWrt, input, sensitivity);
}
