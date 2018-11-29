#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

#include "MUQ/config.h"

#include "MUQ/Modeling/Distributions/Density.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

ImportanceSampling::ImportanceSampling(std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias, pt::ptree const& pt) : SamplingAlgorithm(std::make_shared<SampleCollection>()),
  numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias) {}

ImportanceSampling::ImportanceSampling(std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias, std::vector<Eigen::VectorXd> hyperparameters, pt::ptree const& pt) : SamplingAlgorithm(std::make_shared<SampleCollection>()), numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias), hyperparameters(hyperparameters) {}

std::shared_ptr<SampleCollection> ImportanceSampling::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
  // store a copy of the biasing distribution hyper parameters
  std::vector<Eigen::VectorXd> biasingPara = hyperparameters;

  // insert an empty vector for the state
  biasingPara.insert(biasingPara.begin(), Eigen::VectorXd());

  // loop through the samples to generate them
  for( unsigned int i=0; i<numSamps; ++i ) {
    // sample from the biasing distribution
    biasingPara[0] = bias->Sample(hyperparameters);

    // compute the weight
    const double logweight = target->Evaluate(biasingPara[0])[0](0) - bias->LogDensity(biasingPara);

    // store the sample
    samples->Add(std::make_shared<SamplingState>(biasingPara[0], std::exp(logweight)));
  }

  return samples;
}
