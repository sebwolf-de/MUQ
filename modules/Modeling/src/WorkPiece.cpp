#include "MUQ/Modeling/WorkPiece.h"

// define the muq namespace
using namespace muq::Modeling;

// Create a muq::Modeling::WorkPiece with no fixed number of inputs and outputs and variable input/output types.
WorkPiece::WorkPiece() : 
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(-1), // the number of ouputs is unfixed
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
WorkPiece::WorkPiece(unsigned int const num, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? num : -1), // possibly fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? num : -1), // possibly fix the number of outputs
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
WorkPiece::WorkPiece(unsigned int const numIns, unsigned int const numOuts) : 
  numInputs(numIns), // fix the number of inputs
  numOutputs(numOuts), // fix the number of outputs
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types
WorkPiece::WorkPiece(std::vector<std::string> const& types, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? types.size() : -1), // possibly fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? types.size() : -1), // possibly fix the number of outputs 
  inputTypes(fix==WorkPiece::Fix::Inputs? Types(types) : std::map<unsigned int, std::string>()), // possibly fix the input types
  outputTypes(fix==WorkPiece::Fix::Outputs? Types(types) : std::map<unsigned int, std::string>()), // possibly fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where either some of the inputs have specified types or some of the outputs have specified types
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& types, WorkPiece::Fix const fix) :
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(-1), // the number of ouputs is unfixed
  inputTypes(fix==WorkPiece::Fix::Inputs? types : std::map<unsigned int, std::string>()), // possibly fix the input types
  outputTypes(fix==WorkPiece::Fix::Outputs? types : std::map<unsigned int, std::string>()), // possibly fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types. The number of outputs/inputs (which ever does not have fixed types) is fixed but the types may vary.
WorkPiece::WorkPiece(std::vector<std::string> const& types, unsigned int const num, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? types.size() : num), // fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? types.size() : num), // fix the number of outputs 
  inputTypes(fix==WorkPiece::Fix::Inputs? Types(types) : std::map<unsigned int, std::string>()), // possibly fix the input types
  outputTypes(fix==WorkPiece::Fix::Outputs? Types(types) : std::map<unsigned int, std::string>()), // possibly fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with a fixed number of inputs and outputs with specified types
WorkPiece::WorkPiece(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes) : 
  numInputs(inTypes.size()), // fix the number of inputs 
  numOutputs(outTypes.size()), // fix the number of outputs 
  inputTypes(Types(inTypes)), // fix the input types
  outputTypes(Types(outTypes)), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

std::map<unsigned int, std::string> WorkPiece::Types(std::vector<std::string> const& typesVec) const {
  // initialize the map from input/output number to input/output type
  std::map<unsigned int, std::string> typesMap;

  // populate the map with the elments in the type vector
  for( unsigned int i=0; i<typesVec.size(); ++i ) {
    typesMap[i] = typesVec.at(i);
  }

  return typesMap;
}

unsigned int WorkPiece::SetID() {
  static unsigned int workPieceId = 0;
  return ++workPieceId;
}

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

std::string WorkPiece::Name() const {
  int status;
  std::stringstream ss;

  // the unique name is the name of the (child) class + "_{ID number}"
  ss << abi::__cxa_demangle(typeid(*this).name(), 0, 0, &status) << "_" << id;

  return ss.str();
}

std::string WorkPiece::InputType(unsigned int inputNum) const {
  // make sure the inputNum is less than the number of inputs or that we don't know the number of inputs
  assert(numInputs<0 || inputNum<numInputs);

  if( inputTypes.size()>0 ) { // if we know the input types ...
    int status;
    // ... returned the demangled name
    return abi::__cxa_demangle(inputTypes.at(inputNum).c_str(), 0, 0, &status);
  }

  return "";
}

std::string WorkPiece::OutputType(unsigned int outputNum) const {
  // make sure the inputNum is less than the number of inputs or that we don't know the number of inputs
  assert(numOutputs<0 || outputNum<numOutputs);

  if( outputTypes.size()>0 ) { // if we know the input types ...
    int status;
    // ... returned the demangled name
    return abi::__cxa_demangle(outputTypes.at(outputNum).c_str(), 0, 0, &status);
  }

  return "";
}
