#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

Distribution::~Distribution() {}

Distribution::Distribution() : WorkPiece() {}

std::shared_ptr<Density> Distribution::AsDensity()
{
  return std::make_shared<Density>(shared_from_this());
}

std::shared_ptr<RandomVariable> Distribution::AsVariable()
{
  return std::make_shared<RandomVariable>(shared_from_this());
}

void Distribution::EvaluateImpl(ref_vector<boost::any> const& inputs) {

  // get the mode
  const Distribution::Mode mode = boost::any_cast<Distribution::Mode const>(inputs[0]);

  // the output is either the log-density or a sample
  outputs.resize(1);

  switch( mode ) { // are we evaluting the log density or sampling?
  case Distribution::Mode::EvaluateLogDensity: { // if we are evaluating the log density ...
    // ... store the log density
    outputs[0] = LogDensity(ref_vector<boost::any>(inputs.begin()+1, inputs.end()));
    break;
  }
  case Distribution::Mode::SampleDistribution: { // if we are sampling the distribution
    // .. store the sample
    outputs[0] = Sample(ref_vector<boost::any>(inputs.begin()+1, inputs.end()));
    break;
  }
  default: {
    // something went wrong ...
    assert(false);
  }
  }
}

double Distribution::LogDensity(ref_vector<boost::any> const& inputs) {
  // the first input is always whether we are evaluting the log-density or sampling; we either have n-1 inputs or they are unknown
  assert(numInputs-1==inputs.size() || numInputs<0);

  return LogDensityImpl(inputs);
}

// default behavior of log-density is to return infinity
double Distribution::LogDensityImpl(ref_vector<boost::any> const& inputs) {

  throw muq::NotImplementedError("LogDensityImpl was not implemented for class " + std::string(typeid(*this).name()) );
  return -std::numeric_limits<double>::infinity();
}

std::vector<std::string> Distribution::AddModeInput(std::vector<std::string> const& types) {
  // copy the old types
  std::vector<std::string> new_types(types);

  // the first input is the mode
  new_types.insert(new_types.begin(), typeid(Distribution::Mode).name());

  // return the new types
  return new_types;
}

std::map<unsigned int, std::string> Distribution::AddModeInput(std::map<unsigned int, std::string> const& types) {
  // the first input is the mode
  std::map<unsigned int, std::string> new_types;
  new_types[0] = typeid(Distribution::Mode).name();

  // store all the other types
  for( auto it : types ) {
    new_types[it.first+1] = it.second;
  }

  return new_types;
}

boost::any Distribution::Sample(ref_vector<boost::any> const& inputs) {
  // the first input is always whether we are evaluting the log-density or sampling; we either have n-1 inputs or they are unknown
  assert(numInputs-1==inputs.size() || numInputs<0);

  return SampleImpl(inputs);
}

boost::any Distribution::Sample() {
  // the first input is always whether we are evaluting the log-density or sampling; we either have 1 input or they are unknown
  assert(numInputs==1 || numInputs<0);

  return Sample(ref_vector<boost::any>());
}

boost::any Distribution::SampleImpl(ref_vector<boost::any> const& inputs) {

  throw muq::NotImplementedError("SampleImpl was not implemented for class " + std::string(typeid(*this).name()) );
  return boost::none;
}
