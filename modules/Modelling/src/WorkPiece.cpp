#include "MUQ/Modelling/WorkPiece.h"

using namespace muq::Modelling;

WorkPiece::WorkPiece() : numInputs(-1), numOutputs(-1) {}

WorkPiece::WorkPiece(unsigned int const num, bool const fixInput) : numInputs(fixInput? num : -1), numOutputs(fixInput? -1 : num) {}

WorkPiece::WorkPiece(unsigned int const numIns, unsigned int const numOuts) : numInputs(numIns), numOutputs(numOuts) {}

WorkPiece::WorkPiece(std::vector<std::string> const& types, bool const fixInput) : numInputs(fixInput? types.size() : -1), numOutputs(fixInput? -1 : types.size()), inputTypes(fixInput? types : std::vector<std::string>(0)), outputTypes(fixInput? std::vector<std::string>(0) : types) {}

std::vector<boost::any> WorkPiece::Evaluate() {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || inputs.size()==numInputs);

  outputs.clear();

  // evaluate the WorkPiece
  EvaluateImpl();

  // make sure we have the correct number of outputs
  assert(numOutputs<0 || outputs.size()==numOutputs);

  assert(outputTypes.size()==0 || outputTypes.size()==outputs.size());
  for(unsigned int i=0; i<outputTypes.size(); ++i ) {
    assert(outputTypes[i].compare(outputs[i].type().name())==0);
  }

  return outputs;

}

std::vector<boost::any> WorkPiece::Evaluate(std::vector<boost::any> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  inputs.clear();
  inputs = ins;

  return Evaluate();
}
