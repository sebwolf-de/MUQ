#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"
namespace muq {
  namespace SamplingAlgorithms {

    std::shared_ptr<SamplingState> SubsamplingMIProposal::Sample(std::shared_ptr<SamplingState> currentState) {

      // Consider samples' weights in subsampling, as rejects lead to increased
      // weight on previous samples instead of a new one being added
      //auto sampleCollection = coarseChain->GetSamples();
      for (int i = 0; i < subsampling; i++) {
        sampleWeight++;
        while (sampleWeight >= coarseChain->GetSamples()->samples[sampleID]->weight) {
          if (coarseChain->GetSamples()->samples.size() - 1 == sampleID) {
            coarseChain->Sample();
          } else {
            sampleID++;
            sampleWeight = 0;
          }
        }
      }

      return coarseChain->GetSamples()->samples[sampleID];

      // TODO: This is a far cleaner solution; extremely slow though until iterators for SampleCollections are introduced
      /*sampleID += subsampling;
      *  for (int i = 0; i <= sampleID - coarseChain->GetSamples()->size(); i++) {
      *    coarseChain->Sample();
    }

    return coarseChain->GetSamples()->at(sampleID);*/
    }

  }
}
