#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

using namespace muq::SamplingAlgorithms;

SamplingAlgorithm::SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn) : samples(samplesIn) {}

#if MUQ_HAS_PARCER
SamplingAlgorithm::SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn, std::shared_ptr<parcer::Communicator> comm) : samples(samplesIn), comm(comm) {};
#endif

std::shared_ptr<SampleCollection> SamplingAlgorithm::GetSamples() const { return samples; }

std::shared_ptr<SampleCollection> SamplingAlgorithm::Run() { return Run(std::vector<Eigen::VectorXd>()); }

std::shared_ptr<SampleCollection> SamplingAlgorithm::Run(Eigen::VectorXd const& x0) {return Run(std::vector<Eigen::VectorXd>(1,x0)); }

std::shared_ptr<SampleCollection> SamplingAlgorithm::Run(std::vector<Eigen::VectorXd> const& x0) { return RunImpl(x0); }

#if MUQ_HAS_PARCER
std::shared_ptr<parcer::Communicator> SamplingAlgorithm::GetCommunicator() const { return comm; }
#endif
