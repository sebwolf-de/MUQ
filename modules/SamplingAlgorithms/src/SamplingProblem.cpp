#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingProblem::SamplingProblem(std::shared_ptr<Distribution> target) : target(target) {}

SamplingProblem::~SamplingProblem() {}

boost::any SamplingProblem::SampleTarget(ref_vector<boost::any> const& inputs) const {
  return target->Sample(inputs);
}
