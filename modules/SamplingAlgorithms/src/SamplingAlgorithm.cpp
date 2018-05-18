#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingAlgorithm::SamplingAlgorithm(bool const correlated) : WorkPiece(std::map<unsigned int, std::string>({
      std::pair<unsigned int, std::string>(0, typeid(pt::ptree).name()), // algorithm parameters
	std::pair<unsigned int, std::string>(1, typeid(std::shared_ptr<SamplingProblem>).name()) // sampling problem
	})), correlated(correlated) {}

SamplingAlgorithm::~SamplingAlgorithm() {}

void SamplingAlgorithm::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // get the options for this algorithm
  pt::ptree pt = boost::any_cast<pt::ptree>(inputs.at(0));

  // get the problem that defines the distribution we are trying to characterize/sample
  std::shared_ptr<SamplingProblem> problem = boost::any_cast<std::shared_ptr<SamplingProblem> >(inputs.at(1));
  assert(problem);

  // get the kernel that generates the next sample
  std::shared_ptr<TransitionKernel> kernel = ConstructKernel(pt, problem);
  assert(kernel);

  // get the number of samples
  const unsigned int T = pt.get<unsigned int>("SamplingAlgorithm.NumSamples");

  // the result will be a vector of states
  outputs.resize(1);
  outputs.at(0) = std::vector<std::shared_ptr<SamplingState> >();
  std::vector<std::shared_ptr<SamplingState> >& samples = boost::any_cast< std::vector<std::shared_ptr<SamplingState>>& >(outputs.at(0));
  samples.reserve(T+1);

  // if the samples are correlated, we need a starting sample
  std::shared_ptr<SamplingState> currentState;
  if( correlated ) {
    currentState = std::make_shared<SamplingState>(inputs.at(2), 1.0);
    samples.push_back(currentState);
    kernel->PostStep(1, samples.at(0));
  }

  for(unsigned int t = samples.size(); t<T; ++t) { // loop through each sample

    // get the next sample
    std::shared_ptr<SamplingState> result = kernel->Step(currentState);
  //  const std::vector<boost::any>& result = kernel->Evaluate(kernelInputs);

    // save it if the next sample is valid
    if( currentState != result ) {
       currentState = result;
       samples.push_back(currentState);
     }

     kernel->PostStep(t+1, currentState);

    // allow the kernel to adapt given the step
  }

  // normalize the sample weights
  ReweightSamples(samples);
}

void SamplingAlgorithm::ReweightSamples(std::vector<std::shared_ptr<SamplingState> >& samples) const {
  double totalWeight = 0.0;
  for( unsigned int t=0; t<samples.size(); ++t ) { // loop through each sample
    totalWeight += samples[t]->weight;
  }

  // normalize the weights
  for( unsigned int t=0; t<samples.size(); ++t ) { // loop through each sample
    samples[t]->weight /= totalWeight;
  }
}

boost::any SamplingAlgorithm::FirstMoment() const {
  assert(outputs.size()>0);

  const std::vector<std::shared_ptr<SamplingState> >& samples = boost::any_cast<std::vector<std::shared_ptr<SamplingState> > const&>(outputs[0]);

  return FirstMoment(samples);
}

boost::any SamplingAlgorithm::FirstMoment(std::vector<std::shared_ptr<SamplingState> > const& samples) const {
  // make sure we have an algebra
  assert(algebra);

  // we need at least one sample
  assert(samples.size()>0);

  // compute the weighted sum
  boost::any mean = algebra->Multiply(samples[0]->weight, samples[0]->state);
  for( unsigned int i=1; i<samples.size(); ++i ) {
    mean = algebra->Add(mean, algebra->Multiply(samples[i]->weight, samples[i]->state));
  }

  // return the mean
  return mean;
}
