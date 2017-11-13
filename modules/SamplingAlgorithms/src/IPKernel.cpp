#include "MUQ/SamplingAlgorithms/IPKernel.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(IPKernel)

IPKernel::IPKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem) {}

IPKernel::~IPKernel() {}
