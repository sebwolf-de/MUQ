#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/Utilities/StringUtilities.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

SingleChainMCMC::SingleChainMCMC(boost::property_tree::ptree&             pt,
                                 std::shared_ptr<AbstractSamplingProblem> problem)
{
  numSamps = pt.get<unsigned>("NumSamples");

  std::string kernelString = pt.get<std::string>("KernelList");

  std::vector<std::string> kernelNames = StringUtilities::Split(kernelString, ',');

  unsigned int numBlocks = kernelNames.size();
  kernels.resize(numBlocks);

  // Add the block id to a child tree and construct a kernel for each block
  for(int i=0; i<kernels.size(); ++i){
    boost::property_tree::ptree subTree = pt.get_child(kernelNames.at(i));
    subTree.put("BlockIndex",i);

    kernels.at(i) = TransitionKernel::Construct(subTree, problem);
  }

}

SampleCollection const& SingleChainMCMC::RunImpl(std::vector<boost::any> const& x0)
{
  return samples;
}
