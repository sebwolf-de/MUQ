#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target) : SamplingProblem(target) {}

double ExpensiveSamplingProblem::LogDensity(std::shared_ptr<SamplingState> state) {
  return target->Evaluate(state->state).at(0)(0);
}
