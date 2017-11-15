#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

using namespace muq::SamplingAlgorithms;

ImportanceSampling::ImportanceSampling() :
  SamplingAlgorithm(false) // the sampling algorithm produces uncorrelated samples
{}

ImportanceSampling::~ImportanceSampling() {}

std::shared_ptr<TransitionKernel> ImportanceSampling::Kernel(boost::property_tree::ptree& pt, std::shared_ptr<SamplingProblem> problem) const {
  // we have to use the importance sample transition kernel
  pt.put<std::string>("SamplingAlgorithm.TransitionKernel", "IPKernel");
  
  return TransitionKernel::Construct(pt, problem);
}
