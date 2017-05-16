#include "MUQ/Modeling/WorkPiece.h"

#include <Eigen/Core>

// define the muq namespace
using namespace muq::Modeling;

// Create a muq::Modeling::WorkPiece with no fixed number of inputs and outputs and variable input/output types.
WorkPiece::WorkPiece() : 
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(-1), // the number of ouputs is unfixed
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
WorkPiece::WorkPiece(int const num, WorkPiece::Fix const fix) : 
  numInputs(fix==WorkPiece::Fix::Inputs? num : -1), // possibly fix the number of inputs 
  numOutputs(fix==WorkPiece::Fix::Outputs? num : -1), // possibly fix the number of outputs
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
WorkPiece::WorkPiece(int const numIns, int const numOuts) : 
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
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& types, int const num, WorkPiece::Fix const fixTypes, WorkPiece::Fix const fixNum) :
  numInputs(fixNum==WorkPiece::Fix::Inputs? num : -1), // possibly fix the number of inputs 
  numOutputs(fixNum==WorkPiece::Fix::Outputs? num : -1), // possibly fix the number of outputs 
  inputTypes(fixTypes==WorkPiece::Fix::Inputs? types : std::map<unsigned int, std::string>()), // possibly fix the input types
  outputTypes(fixTypes==WorkPiece::Fix::Outputs? types : std::map<unsigned int, std::string>()), // possibly fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}  

// Create a muq::Modeling::WorkPiece with a fixed number of inputs with specified types and a fixed number of outputs (of uknown type)
WorkPiece::WorkPiece(std::vector<std::string> const& types, int const num) : 
  numInputs(types.size()), // fix the number of inputs 
  numOutputs(num), // fix the number of outputs 
  inputTypes(Types(types)), // fix the input types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece with a fixed number of outputs with specified types and a fixed number of inputs (of uknown type)
WorkPiece::WorkPiece(int const num, std::vector<std::string> const& types) : 
  numInputs(num), // fix the number of inputs 
  numOutputs(types.size()), // fix the number of outputs 
  outputTypes(Types(types)), // fix the input types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the inputs are known and we know the input and output numbers
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const numIns, int const numOuts) :
  numInputs(numIns), // fix the number inputs
  numOutputs(numOuts), // fix the number of outputs
  inputTypes(inTypes), // fix the input types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the outputs are known and we know the input and output numbers
WorkPiece::WorkPiece(int const numIns, std::map<unsigned int, std::string> const& outTypes, int const numOuts) :
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
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const num, std::vector<std::string> const& outTypes) :
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
WorkPiece::WorkPiece(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, int const num) :
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
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const numIn, std::map<unsigned int, std::string> const& outTypes) :
  numInputs(numIn), // fix the number of inputs
  numOutputs(-1), // the number of outputs is unfixed
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of outputs
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, int const numOut) :
  numInputs(-1), // the number of inputs is unfixed
  numOutputs(numOut), // fix the number of outputs
  inputTypes(inTypes), // fix the input types
  outputTypes(outTypes), // fix the output types
  id(SetID()) // the unique id of this WorkPiece
{}

// Create a muq::Modeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of inputs and outputs
WorkPiece::WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const numIn, std::map<unsigned int, std::string> const& outTypes, int const numOut) :
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
  Clear();

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

std::vector<boost::any> WorkPiece::Evaluate(ref_vector<boost::any> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  // we have new outputs
  Clear();

  // the inputs are set, so call evaluate with no inputs
  EvaluateImpl(ins);

  // make sure the output types are correct
  assert(numOutputs<0 || outputs.size()==numOutputs);

  // check the output types
  for( unsigned int i=0; i<outputs.size(); ++i ) {
    assert(CheckOutputType(i, outputs.at(i).type().name()));
  }
  
  return outputs;
}

std::vector<boost::any> WorkPiece::Evaluate(std::vector<boost::any> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  // make sure the input types are correct
  assert(inputTypes.size()==0 || inputTypes.size()==ins.size());
  for(unsigned int i=0; i<inputTypes.size(); ++i ) {
    assert(CheckInputType(i, ins[i].type().name()));
  }
  
  return Evaluate(ToRefVector(ins));
}

boost::any WorkPiece::Jacobian(unsigned int const wrtIn, unsigned int const wrtOut, std::vector<boost::any> const& ins) {
    // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  // make sure the input types are correct
  assert(inputTypes.size()==0 || inputTypes.size()==ins.size());
  for(unsigned int i=0; i<inputTypes.size(); ++i ) {
    assert(CheckInputType(i, ins[i].type().name()));
  }
  
  return Jacobian(wrtIn, wrtOut, ToRefVector(ins));
}

boost::any WorkPiece::Jacobian(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  // clear the outputs and derivative information
  Clear();

  // the inputs are set, so call evaluate with no inputs
  JacobianImpl(wrtIn, wrtOut, ins);

  // make sure the jacobian was computed (the optional jacobian member has a value)
  if( !jacobian ) {
    std::cerr << std::endl << "ERROR: The Jacobian was not computed properly, make sure JacobianImpl gives muq::Modeling::WorkPiece::jacobian a value" << std::endl << std::endl;
    assert(jacobian);
  }

  // return the jacobian (use the * operator because it is a boost::optional)
  return *jacobian;
}

void WorkPiece::JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) {
  // the name of the Eigen::VectorXd's
  const std::string eigenType = typeid(Eigen::VectorXd).name();

  // if both the input and output type is Eigen::VectorXd, default to finite difference
  if( InputType(wrtIn, false).compare(eigenType)==0 && OutputType(wrtOut, false).compare(eigenType)==0 ) {
    std::cout << "USE FINITE DIFFERENCE" << std::endl;
  }

  // invalid! The user has not implemented the Jacobian
  std::cerr << std::endl << "ERROR: No JacobianImpl function for muq::Modeling::WorkPiece implemented, cannot compute Jacobian" << std::endl << std::endl;
  assert(false);
}

boost::any WorkPiece::JacobianAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, std::vector<boost::any> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);

  // clear the outputs and derivative information
  Clear();

  // make sure the input types are correct
  assert(inputTypes.size()==0 || inputTypes.size()==ins.size());
  for(unsigned int i=0; i<inputTypes.size(); ++i ) {
    assert(CheckInputType(i, ins[i].type().name()));
  }
  
  return JacobianAction(wrtIn, wrtOut, vec, ToRefVector(ins));
}

boost::any WorkPiece::JacobianAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& ins) {
  // make sure we have the correct number of inputs
  assert(numInputs<0 || ins.size()==numInputs);
  
  // the inputs are set, so call evaluate with no inputs
  JacobianActionImpl(wrtIn, wrtOut, vec, ins);
  
  // make sure the jacobian was computed (the optional jacobian member has a value)
  if( !jacobianAction ) {
  std::cerr << std::endl << "ERROR: The Jacobian was not computed properly, make sure JacobianActionImpl gives muq::Modeling::WorkPiece::jacobianAction a value" << std::endl << std::endl;
  assert(jacobianAction);
  }

  // return the jacobian action (use the * operator because it is a boost::optional)
  return *jacobianAction;
}

void WorkPiece::JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
   // the name of the Eigen::VectorXd's
  const std::string eigenType = typeid(Eigen::VectorXd).name();

  // if both the input and output type is Eigen::VectorXd, default to finite difference
  if( eigenType.compare(vec.type().name())==0 && InputType(wrtIn, false).compare(eigenType)==0 && OutputType(wrtOut, false).compare(eigenType)==0 ) {
    std::cout << "USE FINITE DIFFERENCE" << std::endl;
  }

  // invalid! The user has not implemented the JacobianAction
  std::cerr << std::endl << "ERROR: No JacobianActionImpl function for muq::Modeling::WorkPiece implemented, cannot compute the action of the Jacobian" << std::endl << std::endl;
  assert(false);
}

std::string WorkPiece::Name() const {
  int status;
  std::stringstream ss;

  // the unique name is the name of the (child) class + "_{ID number}"
  ss << abi::__cxa_demangle(typeid(*this).name(), 0, 0, &status) << "_" << id;

  return ss.str();
}

std::string WorkPiece::InputType(unsigned int inputNum, bool const demangle) const {
  // make sure the inputNum is less than the number of inputs or that we don't know the number of inputs
  assert(numInputs<0 || inputNum<numInputs);

  // an iterator to the input type
  auto it = inputTypes.find(inputNum);

  // we don't know the input type
  if( it==inputTypes.end() ) {
    return "";
  }

  // make it human readable
  if( demangle ) {
    return boost::core::demangle(it->second.c_str());
  }

  // return the input type
  return it->second;
}

std::string WorkPiece::OutputType(unsigned int outputNum, bool const demangle) const {
  // make sure the outputNum is less than the number of outputs or that we don't know the number of outputs
  assert(numOutputs<0 || outputNum<numOutputs);

  // an iterator to the output type
  auto it = outputTypes.find(outputNum);

  // we don't know the output type
  if( it==outputTypes.end() ) {
    return "";
  }

  // make it human readable
  if( demangle ) {
    return boost::core::demangle(it->second.c_str());
  }

  // return the output type
  return it->second;
}

unsigned int WorkPiece::ID() const{
  return id;
}

std::map<unsigned int, std::string> WorkPiece::OutputTypes() const {
  return outputTypes;
}

std::map<unsigned int, std::string> WorkPiece::InputTypes() const {
  return inputTypes;
}

ref_vector<const boost::any> WorkPiece::ToRefVector(std::vector<boost::any> const& anyVec) const {
           
  ref_vector<const boost::any> refs;
  refs.reserve(anyVec.size());

  // populate the input vector
  for(int i=0; i<anyVec.size(); ++i)
    refs.push_back(std::cref(anyVec.at(i)));

  return refs;
}

void WorkPiece::Clear() {
  // we have new outputs --- some WorkPiece's have the outputs stored in a child class so we don't always want to clear them
  if( clearOutputs ) { outputs.clear(); }

  // clear the jacobian
  jacobian = boost::none;

  // clear the jacobian action
  jacobianAction = boost::none;
}

bool WorkPiece::CheckInputType(unsigned int const inputNum, std::string const& type) const {
  // find the input type
  auto it = inputTypes.find(inputNum);
  
  // check to see that the types match (or that we don't know the type)
  if( it!=inputTypes.end() && it->second.compare(type)!=0 ) {
    std::cerr << std::endl << "ERROR: Input types do not match." << std::endl << "\tGiven input: " << boost::core::demangle(type.c_str()) << ", expected " << boost::core::demangle(it->second.c_str()) << std::endl << std::endl;
    return false;
  }
  
  return true;
}

bool WorkPiece::CheckOutputType(unsigned int const outputNum, std::string const& type) const {
  // find the output type
  auto it = outputTypes.find(outputNum);
  
  // check to see that the types match (or that we don't know the type)
  if( it!=outputTypes.end() && it->second.compare(type)!=0 ) {
    std::cerr << std::endl << "ERROR: Output types do not match." << std::endl << "\tGiven input: " << boost::core::demangle(type.c_str()) << ", expected " << boost::core::demangle(it->second.c_str()) << std::endl << std::endl;
    return false;
  }
  
  return true;
}
