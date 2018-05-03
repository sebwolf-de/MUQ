#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> targetIn,
                                 std::vector<int> const& inputSizes) : AbstractSamplingProblem(inputSizes.size(), inputSizes),
                                                                       target(targetIn) {}

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::Distribution> targetIn) : SamplingProblem(targetIn, GetBlockSizes(targetIn)) {}


double SamplingProblem::LogDensity(std::shared_ptr<SamplingState> state) {
  assert(target);
  return target->LogDensity(state->state);
}

boost::any SamplingProblem::GradLogDensity(std::shared_ptr<SamplingState> state,
                                           unsigned                       blockWrt)
{
  return target->Jacobian(blockWrt, 0, state->state);
}

std::vector<int> SamplingProblem::GetBlockSizes(std::shared_ptr<Distribution> target)
{
  int numBlocks = GetNumBlocks(target);

  std::vector<int> output(numBlocks);
  for(int i=0; i<numBlocks; ++i)
    output.at(i) = target->InputSize(i);

  return output;
}

unsigned SamplingProblem::GetNumBlocks(std::shared_ptr<Distribution> target)
{
  if(target->numInputs < 0){
    throw std::invalid_argument("When not manually specified, \"SamplingProblem\" requires the target distribution to have a specified number of inputs, but \"" + target->Name() + "\" does not specify a fixed number of inputs, i.e., target->numInputs < 0.");
  }
  return target->numInputs;
}
