#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingProblem::SamplingProblem(std::shared_ptr<Distribution> target) : target(target) {}

SamplingProblem::SamplingProblem(std::shared_ptr<Distribution> target, std::shared_ptr<Distribution> bias) : target(target), bias(bias) {}

SamplingProblem::~SamplingProblem() {}

boost::any SamplingProblem::SampleTarget(ref_vector<boost::any> const& inputs) const {
  assert(target);
  return target->Sample(inputs);
}

double SamplingProblem::EvaluateLogTarget(muq::Modeling::ref_vector<boost::any> const& inputs) const {
  assert(target);
  return target->LogDensity(inputs);
}

boost::any SamplingProblem::SampleBiasingDistribution(muq::Modeling::ref_vector<boost::any> const& inputs) const {
  assert(bias);
  return (*bias)->Sample(inputs);
}

double SamplingProblem::EvaluateLogBiasingDistribution(muq::Modeling::ref_vector<boost::any> const& inputs) const {
  assert(bias);
  return (*bias)->LogDensity(inputs);
}
