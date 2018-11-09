#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

#include "MUQ/Modeling/Distributions/Density.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

ImportanceSampling::ImportanceSampling(std::shared_ptr<Density> const& target, std::shared_ptr<Distribution> const& bias, pt::ptree const& pt) : SamplingAlgorithm(std::make_shared<SampleCollection>()),
  numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias) {}

ImportanceSampling::~ImportanceSampling() {}

std::shared_ptr<SampleCollection> ImportanceSampling::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
  // loop through the samples to generate them
  for( unsigned int i=0; i<numSamps; ++i ) {
    const Eigen::VectorXd& state = bias->Sample();
    const double logweight = target->LogDensity(std::vector<Eigen::VectorXd>(1, state)) - bias->LogDensity(state);

    samples->Add(std::make_shared<SamplingState>(state, std::exp(logweight)));
  }

  return samples;
}
