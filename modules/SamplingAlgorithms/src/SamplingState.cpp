#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::SamplingAlgorithms;

SamplingState::SamplingState(boost::any const& state, double const weight) : state(state), weight(weight) {}

SamplingState::~SamplingState() {}

bool SamplingState::HasMeta(std::string const& metaKey)
{
  auto iter = meta.find(metaKey);
  return iter != meta.end();
}
