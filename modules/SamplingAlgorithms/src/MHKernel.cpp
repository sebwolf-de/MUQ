#include "MUQ/SamplingAlgorithms/MHKernel.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MHKernel)

MHKernel::MHKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem) {}

MHKernel::~MHKernel() {}
