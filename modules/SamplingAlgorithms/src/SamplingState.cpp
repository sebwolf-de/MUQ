#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::SamplingAlgorithms;

SamplingState::SamplingState(boost::any const& state, double const weight) : state(state), weight(weight) {}

SamplingState::~SamplingState() {}
