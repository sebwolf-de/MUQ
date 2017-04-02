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

// Create a muq::Modeling::WorkPiece where either some of the inputs have specified types or some of the outputs have specified types and either the number of inputs or the number of outputs is fixed
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& types, unsigned int const num, WorkPiece::Fix const fixTypes, WorkPiece::Fix const fixNum) :
  numInputs(fixNum==WorkPiece::Fix::Inputs? num : -1), // possibly fix the number of inputs 
  numOutputs(fixNum==WorkPiece::Fix::Outputs? num : -1), // possibly fix the number of outputs 
  inputTypes(fixTypes==WorkPiece::Fix::Inputs? types : std::map<unsigned int, std::string>()), // possibly fix the input types
  outputTypes(fixTypes==WorkPiece::Fix::Outputs? types : std::map<unsigned int, std::string>()), // possibly fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}  

// Create a muq::Modeling::WorkPiece with a fixed number of inputs with specified types and a fixed number of outputs (of uknown type)
WorkPiece::WorkPiece(std::vector<std::string> const& types, unsigned int const num) : 
  numInputs(types.size()), // fix the number of inputs 
  numOutputs(num), // fix the number of outputs 
  inputTypes(Types(types)), // fix the input types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with a fixed number of outputs with specified types and a fixed number of inputs (of uknown type)
WorkPiece::WorkPiece(unsigned int const num, std::vector<std::string> const& types) : 
  numInputs(num), // fix the number of inputs 
  numOutputs(types.size()), // fix the number of outputs 
  outputTypes(Types(types)), // fix the input types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the inputs are known and we know the input and output numbers
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, unsigned int const numIns, unsigned int const numOuts) :
  numInputs(numIns), // fix the number inputs
  numOutputs(numOuts), // fix the number of outputs
  inputTypes(inTypes), // fix the input types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the outputs are known and we know the input and output numbers
WorkPiece::WorkPiece(unsigned int const numIns, std::map<unsigned int, std::string> const& outTypes, unsigned int const numOuts) :
  numInputs(numIns), // fix the number inputs
  numOutputs(numOuts), // fix the number of outputs
  outputTypes(outTypes), // fix the output types
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

// Create a muq::Mdoeling::WorkPiece where some of the inputs and all of the outputs have specified types
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::vector<std::string> const& outTypes) :
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(outTypes.size()), // fix the number of outputs
  inputTypes(inTypes), // fix the inputs types
  outputTypes(Types(outTypes)), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the inputs are known with a known number of inputs and all of the outputs have specified types
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, unsigned int const num, std::vector<std::string> const& outTypes) :
  numInputs(num), // fix the number of inputs
  numOutputs(outTypes.size()), // fix the number of outputs
  inputTypes(inTypes), // fix the inputs types
  outputTypes(Types(outTypes)), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Mdoeling::WorkPiece where some of the outputs and all of the inputs have specified types
WorkPiece::WorkPiece(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes) :
  numInputs(inTypes.size()), // fix the number of inputs
  numOutputs(-1), // the number of outputs is unfixed
  inputTypes(Types(inTypes)), // fix the inputs types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the outputs with a known number of outputs and all of the inputs have specified types
WorkPiece::WorkPiece(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, unsigned int const num) :
  numInputs(inTypes.size()), // fix the number of inputs
  numOutputs(num), // fix the number of outputs
  inputTypes(Types(inTypes)), // fix the inputs types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes) :
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(-1), // the number of outputs is unfixed
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of inputs
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, unsigned int const numIn, std::map<unsigned int, std::string> const& outTypes) :
  numInputs(numIn), // fix the number of inputs
  numOutputs(-1), // the number of outputs is unfixed
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of outputs
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, unsigned int const numOut) :
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(numOut), // fix the number of outputs
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of inputs and outputs
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, unsigned int const numIn, std::map<unsigned int, std::string> const& outTypes, unsigned int const numOut) :
  numInputs(numIn), // fix the number of inputs
  numOutputs(numOut), // fix the number of outputs
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes), // fix the output types
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

std::vector<std::string> WorkPiece::Types(std::vector<boost::any> const& vec) const {
    // create a vector of the types
  std::vector<std::string> types;
  types.reserve(vec.size());

  // populate types with the type of each element of vec
  for( auto it : vec ) {
    types.push_back(it.type().name());
  }

  // the types and vector should be the same size
  assert(types.size()==vec.size());

  return types;
}

unsigned int WorkPiece::SetID() {
  static unsigned int workPieceId = 0;
  return ++workPieceId;
}

std::vector<boost::any> WorkPiece::Evaluate() {
  // make sure we have the correct number of inputs
  assert(numInputs<=0);

  // clear the outputs
  if( clearOutputs ) { outputs.clear(); }

  // evaluate the WorkPiece
  std::vector<std::reference_wrapper<const boost::any>> emptyVec;
  EvaluateImpl(emptyVec);

  // make sure we have the correct number of outputs
  assert(numOutputs<0 || outputs.size()==numOutputs);

  // make sure the output types are correct
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
  
  for( unsigned int i=0; i<outputs.size(); ++i ) {
    // find the output type
    auto it = outputTypes.find(i);

    if( it!=outputTypes.end() ) { // if we know the output type
      // check to see that the types match
      assert(it->second.compare(outputs.at(i).type().name())==0);
    }
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
