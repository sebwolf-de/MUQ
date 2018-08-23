#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace pt = boost::property_tree;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, pt::ptree const& pt) : SamplingProblem(target) {
  const std::string& regOpt = pt.get<std::string>("RegressionOptions");
  //reg = std::make_shared<LocalRegression>
}

double ExpensiveSamplingProblem::LogDensity(std::shared_ptr<SamplingState> state) {
  return target->Evaluate(state->state).at(0)(0);
}
