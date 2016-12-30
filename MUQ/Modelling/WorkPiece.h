#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<cassert>
#include<memory>

#include "boost/any.hpp"

namespace muq {
  namespace Modelling {
    /// Base class for MUQ's modelling envronment
    class WorkPiece {
    protected:
      
      /// Does the constructor fix the inputs or the outputs?
      enum Fix { 
	/// The constructor fixes the input number and possibly the types
	Inputs, 
	/// The constructor fixes the output number and possibly the types
	Outputs };
      
    public:
      
      /// Create a muq::Modelling::WorkPiece with no fixed number of inputs and outputs and variable input/output types.
      WorkPiece();

      /// Create a muq::Modelling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
      /**
	 @param[in] num The number of inputs or outputs (which one depends on the second parameter)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the first parameter is the number of inputs; WorkPiece::Fix::Outputs: the first parameter is the number of outputs
       */
      WorkPiece(unsigned int const num, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);
            
      /// Create a muq::Modelling::WorkPiece with a fixed number of inputs and outputs but variable input/output types
      /**
	 @param[in] numIns The number of inputs
	 @param[in] numOuts The number of outputs
       */
      WorkPiece(unsigned int const numIns, unsigned int const numOuts);

      /// Create a muq::Modelling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types
      /**
	 If the number and type of the inputs is specified then the number and type of the outputs is variable.  The opposite is true if the number and type of the outputs is specified.
	 @param[in] types A vector of strings, each element is the type of an input or output (the number of inputs or outputs is the size of this vector)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs
       */
      WorkPiece(std::vector<std::string> const& types, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);

      /// Create a muq::Modelling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types. The number of outputs/inputs (which ever does not have fixed types) is fixed but the types may vary.
      /**
	 If the number and type of the inputs is specified then the number of outputs is fixed but the type of the outputs is variable.  The opposite is true if the number and type of the outputs is specified.
	 @param[in] types A vector of strings, each element is the type of an input or output (the number of inputs or outputs is the size of this vector)
	 @param[in] num The number of outputs or inputs (which one depends on the third parameter)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs and the second parameter is the number of outputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs an the second parameter is the number of inputs
       */
      WorkPiece(std::vector<std::string> const& types, unsigned int const num, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);

      /// Create a muq::Modelling::WorkPiece with a fixed number of inputs and outputs with specified types
      /**
	 @param[in] inTypes A vector of strings, each element is the type of an input (the number of inputs is the size of this vector)
	 @param[in] outTypes A vector of strings, each element is the type of an output (the number of outputs is the size of this vector)
       */
      WorkPiece(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes);
            
      virtual ~WorkPiece() {}

      std::vector<boost::any> Evaluate(std::vector<boost::any> const& ins);

      std::vector<boost::any> Evaluate();
      
      template<typename... Args>			
	std::vector<boost::any> Evaluate(Args... args) {
	inputs.clear();
	outputs.clear();

	inputs.resize((numInputs<0? 0 : numInputs));

	// begin calling the EvaluateMulti with the first input
	return EvaluateMulti(0, args...);
      }	

      std::vector<boost::any> EvaluateMulti(unsigned int const inputNum) {
	return Evaluate();
      }

      template<typename ith, typename... Args>		       
	std::vector<boost::any> EvaluateMulti(unsigned int const inputNum, ith const& in, Args... args) {
	// we have not yet put all of the inputs into the map, the ith should be less than the total number
	assert(numInputs<0 || inputNum+1<numInputs);
	
	if( inputTypes.size()>0 ) { // if we know the input types
	  // make sure the index is valid
	  assert(inputNum+1<inputTypes.size());

	  // make sure the type is correct
	  assert(inputTypes[inputNum].compare(typeid(in).name())==0);
	}

	// add the last input to the input vector
	if( numInputs<0 ) { inputs.push_back(in); } else { inputs[inputNum] = in; }

	// call with EvaluateMulti with the remaining inputs
	return EvaluateMulti(inputNum+1, args...);
      }								
      
      template<typename last>			
	std::vector<boost::any> EvaluateMulti(unsigned int const inputNum, last const& in) {
	// this is the last input, the last one should equal the total number of inputs
	assert(numInputs<0 || inputNum+1==numInputs);

	if( inputTypes.size()>0 ) { // if we know the input types
	  // make sure the index is valid
	  assert(inputNum+1==inputTypes.size());

	  // make sure the type is correct
	  assert(inputTypes[inputNum].compare(typeid(in).name())==0);
	}

	// add the last input to the input vector
	if( numInputs<0 ) { inputs.push_back(in); } else { inputs[inputNum] = in; }

	return Evaluate();
      }
      
      /// The number of inputs
      const int numInputs;

      /// The number of outputs
      const int numOutputs;

    protected:

      std::vector<boost::any> inputs = std::vector<boost::any>(0);

      std::vector<boost::any> outputs = std::vector<boost::any>(0);

    private:

      virtual void EvaluateImpl() = 0;

      std::vector<std::string> inputTypes = std::vector<std::string>(0);

      std::vector<std::string> outputTypes = std::vector<std::string>(0);
      
    }; // class WorkPiece
  } // namespace Modelling
} // namespace muq
