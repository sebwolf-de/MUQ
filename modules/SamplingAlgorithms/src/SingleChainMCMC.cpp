#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

#include "MUQ/SamplingAlgorithms/MarkovChain.h"

#include "MUQ/Utilities/StringUtilities.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

SingleChainMCMC::SingleChainMCMC(boost::property_tree::ptree&             pt,
                                 std::shared_ptr<AbstractSamplingProblem> problem) : SamplingAlgorithm(std::make_shared<MarkovChain>()),
                                                                                     printLevel(pt.get("PrintLevel",3))
{
  numSamps = pt.get<unsigned int>("NumSamples");
  burnIn = pt.get("BurnIn",0);

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

void SingleChainMCMC::PrintStatus(std::string prefix, unsigned int currInd) const
{
  std::cout << prefix << int(std::floor(double((currInd - 1) * 100) / double(numSamps))) << "% Complete" << std::endl;
  if(printLevel>1){
    for(int blockInd=0; blockInd<kernels.size(); ++blockInd){
      std::cout << prefix << "  Block " << blockInd << ":\n";
      kernels.at(blockInd)->PrintStatus(prefix + "    ");
    }
  }
}

std::shared_ptr<SampleCollection> SingleChainMCMC::RunImpl(std::vector<Eigen::VectorXd> const& x0)
{
  unsigned sampNum = 0;
  std::vector<std::shared_ptr<SamplingState>> newStates;
  std::shared_ptr<SamplingState> prevState = std::make_shared<SamplingState>(x0);

  samples->Add(prevState);

  // What is the next iteration that we want to print at
  const unsigned int printIncr = std::floor(numSamps / double(10));
  unsigned int nextPrintInd = printIncr;

  // Run until we've received the desired number of samples
  if(printLevel>0)
    std::cout << "Starting single chain MCMC sampler..." << std::endl;

  auto startTime = std::chrono::high_resolution_clock::now();
  while(sampNum < numSamps)
  {
    // Should we print
    if(sampNum > nextPrintInd){
      if(printLevel>0){
        PrintStatus("  ", sampNum);
      }
      nextPrintInd += printIncr;
    }

    // Loop through each parameter block
    for(int blockInd=0; blockInd<kernels.size(); ++blockInd){

      kernels.at(blockInd)->PreStep(sampNum, prevState);

      newStates = kernels.at(blockInd)->Step(sampNum, prevState);

      // Add the new states to the SampleCollection
      for(int i=0; i<newStates.size(); ++i){
        sampNum++;

        if(newStates.at(i)!=prevState){
          prevState = newStates.at(i);

          if(sampNum>=burnIn)
            samples->Add(newStates.at(i));

        }else{
          prevState->weight += 1;
          if(sampNum==burnIn)
            samples->Add(prevState);
        }
      }

      kernels.at(blockInd)->PostStep(sampNum, newStates);
    }
  }

  auto endTime = std::chrono::high_resolution_clock::now();
  double runTime = std::chrono::duration<double>(endTime - startTime).count();

  if(printLevel>0){
    PrintStatus("  ", numSamps+1);
    std::cout << "Completed in " << runTime << " seconds." << std::endl;
  }


  return samples;
}
