#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TransitionKernel::TransitionKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) :
  WorkPiece(/*std::vector<std::string>(1,typeid(std::shared_ptr<SamplingState>).name()), WorkPiece::Fix::Outputs*/), problem(problem) {}

TransitionKernel::~TransitionKernel() {}

std::shared_ptr<TransitionKernel::TransitionKernelMap> TransitionKernel::GetTransitionKernelMap() {
  // define a static map from type to constructor
  static std::shared_ptr<TransitionKernelMap> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<TransitionKernelMap>();
  }

  return map;
}

std::shared_ptr<TransitionKernel> TransitionKernel::Construct(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) {
  // get the name of the kernel
  const std::string& kernelName = pt.get<std::string>("SamplingAlgorithm.TransitionKernel");

  // construct it from the map
  auto kernelMap = GetTransitionKernelMap();
  auto iter = kernelMap->find(kernelName);
  if(iter == kernelMap->end()){
    std::cerr << "ERROR: Could not find Transition Kernel \"" << kernelName << "\", available types are:\n";
    
    for(auto it=kernelMap->begin(); it!=kernelMap->end(); ++it)
      std::cerr << "  " << it->first << std::endl;
    std::cerr << std::endl;

    assert(iter != kernelMap->end());

  }

  return iter->second(pt, problem);

}

void TransitionKernel::PostStep(unsigned int const t, std::shared_ptr<SamplingState> state) {}
