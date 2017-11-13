#include "MUQ/SamplingAlgorithms/MCKernel.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MCKernel)

MCKernel::MCKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem) {}

MCKernel::~MCKernel() {}
