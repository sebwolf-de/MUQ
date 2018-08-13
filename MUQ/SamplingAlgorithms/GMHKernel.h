#ifndef GMHKERNEL_H_
#define GMHKERNEL_H_

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
    class GMHKernel : public MHKernel {
    public:

      GMHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem);
      
      GMHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn);

      virtual ~GMHKernel();

      virtual void PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) override;

      virtual std::vector<std::shared_ptr<SamplingState> > Step(unsigned int const t, std::shared_ptr<SamplingState> state) override;

    private:
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
