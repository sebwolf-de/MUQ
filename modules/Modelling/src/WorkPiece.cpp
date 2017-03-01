#include "MUQ/Modelling/WorkPiece.h"

// define the muq namespace
using namespace muq::Modelling;

// Create a muq::Modeling::WorkPiece with no fixed number of inputs and outputs and variable input/output types.
WorkPiece::WorkPiece() : 
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(-1) // the number of ouputs is unfixed
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
WorkPiece::WorkPiece(unsigned int const num, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? num : -1), // possibly fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? num : -1) // possibly fix the number of outputs
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
WorkPiece::WorkPiece(unsigned int const numIns, unsigned int const numOuts) : 
  numInputs(numIns), // fix the number of inputs
  numOutputs(numOuts) // fix the number of outputs
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types
WorkPiece::WorkPiece(std::vector<std::string> const& types, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? types.size() : -1), // possibly fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? types.size() : -1), // possibly fix the number of outputs 
  inputTypes(fix==WorkPiece::Fix::Inputs? types : std::vector<std::string>(0)), // possibly fix the input types
  outputTypes(fix==WorkPiece::Fix::Outputs? types : std::vector<std::string>(0)) // possibly fix the output types
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types. The number of outputs/inputs (which ever does not have fixed types) is fixed but the types may vary.
WorkPiece::WorkPiece(std::vector<std::string> const& types, unsigned int const num, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? types.size() : num), // fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? types.size() : num), // fix the number of outputs 
  inputTypes(fix==WorkPiece::Fix::Inputs? types : std::vector<std::string>(0)), // possibly fix the input types
  outputTypes(fix==WorkPiece::Fix::Outputs? types : std::vector<std::string>(0)) // possibly fix the output types
{}

// Create a muq::Modeling::WorkPiece with a fixed number of inputs and outputs with specified types
WorkPiece::WorkPiece(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes) : 
  numInputs(inTypes.size()), // fix the number of inputs 
  numOutputs(outTypes.size()), // fix the number of outputs 
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes) // fix the output types
{}

std::vector<boost::any> WorkPiece::Evaluate() {
  // make sure we have the correct number of inputs
  assert(numInputs<=0);

  outputs.clear();

  // evaluate the WorkPiece
  std::vector<std::reference_wrapper<const boost::any>> emptyVec;
  EvaluateImpl(emptyVec);

  // make sure we have the correct number of outputs
  assert(numOutputs<0 || outputs.size()==numOutputs);

  // make sure the output types are correct
  assert(outputTypes.size()==0 || outputTypes.size()==outputs.size());
  for(unsigned int i=0; i<outputTypes.size(); ++i ) {
    assert(outputTypes[i].compare(outputs[i].type().name())==0);
  }

  // return the outputs
  return outputs;
}

std::vector<boost::any> WorkPiece::Evaluate(std::vector<std::reference_wrapper<const boost::any>> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  // the inputs are set, so call evaluate with no inputs
  EvaluateImpl(ins);

  // make sure the output types are correct
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

  // make sure the input types are correct
  assert(inputTypes.size()==0 || inputTypes.size()==ins.size());
  for(unsigned int i=0; i<inputTypes.size(); ++i ) {
    assert(inputTypes[i].compare(ins[i].type().name())==0);
  }
  
  // create the input vector and reserve enough space
  ref_vector<const boost::any> in_refs;
  in_refs.reserve(ins.size());
  // populate the input vector
  for(int i=0; i<ins.size(); ++i)
    in_refs.push_back(std::cref(ins.at(i)));
  
  // the inputs are set, so call evaluate with no inputs
  return Evaluate(in_refs);
}
