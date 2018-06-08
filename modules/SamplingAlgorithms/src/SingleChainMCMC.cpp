#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/Utilities/StringUtilities.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

SingleChainMCMC::SingleChainMCMC(boost::property_tree::ptree&             pt,
                                 std::shared_ptr<AbstractSamplingProblem> problem)
{
  numSamps = pt.get<unsigned int>("NumSamples");
  burnIn = pt.get("BurnIn",0);

  std::string kernelString = pt.get<std::string>("KernelList");

  std::vector<std::string> kernelNames = StringUtilities::Split(kernelString, ',');

  unsigned int numBlocks = kernelNames.size();
  kernels.resize(numBlocks);

  scheduler = std::make_shared<ThinScheduler>(pt);

  // Add the block id to a child tree and construct a kernel for each block
  for(int i=0; i<kernels.size(); ++i){
    boost::property_tree::ptree subTree = pt.get_child(kernelNames.at(i));
    subTree.put("BlockIndex",i);

    kernels.at(i) = TransitionKernel::Construct(subTree, problem);
  }

}

SampleCollection const& SingleChainMCMC::RunImpl(std::vector<Eigen::VectorXd> const& x0)
{
  unsigned sampNum = 0;
  std::vector<std::shared_ptr<SamplingState>> newStates;
  std::shared_ptr<SamplingState> prevState = std::make_shared<SamplingState>(x0);

  samples.Add(prevState);

  // Run until we've received the desired number of samples
  while(sampNum < numSamps)
  {
    // Loop through each parameter block
    for(int blockInd=0; blockInd<kernels.size(); ++blockInd){

      kernels.at(blockInd)->PreStep(sampNum, prevState);

      newStates = kernels.at(blockInd)->Step(sampNum, prevState);

      // Add the new states to the SampleCollection
      for(int i=0; i<newStates.size(); ++i){
        sampNum++;

        if(newStates.at(i)!=prevState){
          prevState = newStates.at(i);

          if((sampNum>=burnIn)&&(scheduler->ShouldSave(sampNum)))
            samples.Add(newStates.at(i));

        }else{
          prevState->weight += 1;
          if((sampNum==burnIn)&&(scheduler->ShouldSave(sampNum)))
            samples.Add(prevState);
        }
      }

      kernels.at(blockInd)->PostStep(sampNum, newStates);
    }
  }

  return samples;
}
