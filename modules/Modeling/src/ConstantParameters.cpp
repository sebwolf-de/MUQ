#include "MUQ/Modeling/ConstantParameters.h"

using namespace muq::Modeling;

// create a WorkPiece with no inputs, known output types, and known output number
ConstantParameters::ConstantParameters(std::vector<boost::any> const& outs) : WorkPiece(0, Types(outs)), outs(outs) {
  // the outputs will never change so we should not clear them
  clearOutputs = false;

  // populate the vector of outputs
  outputs.resize(numOutputs);
  std::copy(outs.begin(), outs.end(), outputs.begin());
}

// the outputs are already set and not cleared so don't do anything
void ConstantParameters::EvaluateImpl(ref_vector<boost::any> const& inputs) {}



