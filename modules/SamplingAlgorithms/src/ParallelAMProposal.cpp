#include "MUQ/SamplingAlgorithms/ParallelAMProposal.h"

#if MUQ_HAS_PARCER

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_MCMC_PROPOSAL(ParallelAMProposal)
ParallelAMProposal::ParallelAMProposal(pt::ptree const& pt , std::shared_ptr<AbstractSamplingProblem> prob) : AMProposal(pt, prob) {}

#endif // end MUQ_HAS_PARCER
