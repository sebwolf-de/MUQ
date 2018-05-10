#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::WorkPiece> targetIn,
                                 std::vector<int> const& inputSizes) : AbstractSamplingProblem(inputSizes.size(), inputSizes),
                                                                       target(targetIn) {}

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::WorkPiece> targetIn) : SamplingProblem(targetIn, GetBlockSizes(targetIn)) {}


double SamplingProblem::LogDensity(std::shared_ptr<SamplingState> state) {
  assert(target);

  std::vector<boost::any> anyState(state->state.size());
  for(int i=0; i<state->state.size(); ++i)
    anyState.at(i) = state->state.at(i);

  return boost::any_cast<double>(target->Evaluate(anyState).at(0));
}

Eigen::VectorXd SamplingProblem::GradLogDensity(std::shared_ptr<SamplingState> state,
                                           unsigned                       blockWrt)
{
  return boost::any_cast<Eigen::MatrixXd>(target->Jacobian(blockWrt, 0, state->state));
}

std::vector<int> SamplingProblem::GetBlockSizes(std::shared_ptr<WorkPiece> target)
{
  int numBlocks = GetNumBlocks(target);

  std::vector<int> output(numBlocks);
  for(int i=0; i<numBlocks; ++i)
    output.at(i) = target->InputSize(i);

  return output;
}

unsigned SamplingProblem::GetNumBlocks(std::shared_ptr<WorkPiece> target)
{
  if(target->numInputs < 0){
    throw std::invalid_argument("When not manually specified, \"SamplingProblem\" requires the target distribution to have a specified number of inputs, but \"" + target->Name() + "\" does not specify a fixed number of inputs, i.e., target->numInputs < 0.");
  }
  return target->numInputs;
}
