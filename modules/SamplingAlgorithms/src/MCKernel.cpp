#include "MUQ/SamplingAlgorithms/MCKernel.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MCKernel)

MCKernel::MCKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem), N(pt.get<unsigned int>("SamplingAlgorithm.NumSamples")) {}

MCKernel::~MCKernel() {}

void MCKernel::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  const boost::any state = problem->SampleTarget(inputs);

  outputs.resize(1);
  outputs[0] = std::make_shared<SamplingState>(state, 1.0);
}
