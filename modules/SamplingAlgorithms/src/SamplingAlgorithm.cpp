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
  pt::ptree pt = boost::any_cast<pt::ptree>(inputs[0]);

  // get the problem that defines the distribution we are trying to characterize/sample
  std::shared_ptr<SamplingProblem> problem = boost::any_cast<std::shared_ptr<SamplingProblem> >(inputs[1]);
  assert(problem);

  // get the kernel that generates the next sample
  std::shared_ptr<TransitionKernel> kernel = Kernel(pt, problem);

  // get the number of samples
  const unsigned int T = pt.get<unsigned int>("SamplingAlgorithm.NumSamples");

  // the result will be a vector of states
  outputs.resize(1);
  outputs[0] = std::vector<std::shared_ptr<SamplingState> >();
  std::vector<std::shared_ptr<SamplingState> >& samples = boost::any_cast<std::vector<std::shared_ptr<SamplingState> >&>(outputs[0]);
  samples.reserve(T+1);

  // if the samples are correlated, we need a starting sample
  boost::any init;
  if( correlated ) {
    init = std::make_shared<SamplingState>(inputs[2], 1.0);
    samples.push_back(boost::any_cast<std::shared_ptr<SamplingState> >(init));

    // begin the kernel adaptation
    kernel->PostStep(1, samples[0]);
  }

  // copy the inputs---these are the inputs to the distributions (e.g., target, biasing, and/or propsal)
  ref_vector<boost::any> distInputs(inputs.begin()+(correlated?3:2), inputs.end());
  if( correlated ) { distInputs.insert(distInputs.begin(), std::cref(init)); }

  unsigned int t = samples.size();
  while( t<T ) { // loop through each sample
    // get the next sample
    const std::vector<boost::any>& result = kernel->Evaluate(distInputs);

    // save it if the next sample is valid
    if( result[0].type()==typeid(std::shared_ptr<SamplingState>) ) { samples.push_back(boost::any_cast<std::shared_ptr<SamplingState> >(result[0])); }

    // allow the kernel to adapt given the step
    kernel->PostStep(++t, samples[samples.size()-1]);

    // Update the current sample in the correlated case
    if( correlated ) {
      init = samples[samples.size()-1];
      distInputs[0] = std::cref(init);
    }
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
