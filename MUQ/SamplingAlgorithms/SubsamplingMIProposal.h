#ifndef SUBSAMPLINGMIPROPOSAL_H
#define SUBSAMPLINGMIPROPOSAL_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {
    class SubsamplingMIProposal : public MCMCProposal {
    public:
      SubsamplingMIProposal (pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<SingleChainMCMC> coarseChain)
      : MCMCProposal(pt,prob), coarseChain(coarseChain),
        subsampling(pt.get("subsampling",1))
      {}

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) override {
        return 0;
      }

    private:
      std::shared_ptr<SingleChainMCMC> coarseChain;
      int sampleID = 0;
      int sampleWeight = 0;
      const int subsampling;
    };
  }
}

#endif
