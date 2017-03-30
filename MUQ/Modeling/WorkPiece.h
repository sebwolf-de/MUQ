#ifndef WORKPIECE_H_
#define WORKPIECE_H_

#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<cassert>
#include<memory>

#include "boost/any.hpp"

namespace muq {
  namespace Modeling {
    /// A vector of references to something ... 
    template <typename T>
      using ref_vector = std::vector<std::reference_wrapper<const T>>;
    
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
      
      /// Create a muq::Modeling::WorkPiece with no fixed number of inputs and outputs and variable input/output types.
      WorkPiece();
      
      /// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
      /**
	 @param[in] num The number of inputs or outputs (which one depends on the second parameter)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the first parameter is the number of inputs; WorkPiece::Fix::Outputs: the first parameter is the number of outputs
      */
      WorkPiece(unsigned int const num, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);
      
      /// Create a muq::Modeling::WorkPiece with a fixed number of inputs and outputs but variable input/output types
      /**
	 @param[in] numIns The number of inputs
	 @param[in] numOuts The number of outputs
      */
      WorkPiece(unsigned int const numIns, unsigned int const numOuts);
      
      /// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types
      /**
	 If the number and type of the inputs is specified then the number and type of the outputs is variable.  The opposite is true if the number and type of the outputs is specified.
	 @param[in] types A vector of strings, each element is the type of an input or output (the number of inputs or outputs is the size of this vector)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs
      */
      WorkPiece(std::vector<std::string> const& types, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);
      
      /// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types. The number of outputs/inputs (which ever does not have fixed types) is fixed but the types may vary.
      /**
	 If the number and type of the inputs is specified then the number of outputs is fixed but the type of the outputs is variable.  The opposite is true if the number and type of the outputs is specified.
	 @param[in] types A vector of strings, each element is the type of an input or output (the number of inputs or outputs is the size of this vector)
	 @param[in] num The number of outputs or inputs (which one depends on the third parameter)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs and the second parameter is the number of outputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs an the second parameter is the number of inputs
      */
      WorkPiece(std::vector<std::string> const& types, unsigned int const num, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);
      
      /// Create a muq::Modeling::WorkPiece with a fixed number of inputs and outputs with specified types
      /**
	 @param[in] inTypes A vector of strings, each element is the type of an input (the number of inputs is the size of this vector)
	 @param[in] outTypes A vector of strings, each element is the type of an output (the number of outputs is the size of this vector)
      */
      WorkPiece(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes);
      
      /// Default destructor
      virtual ~WorkPiece() {}
      
      /// Evaluate this muq::Modeling::WorkPiece
      /**
	 This function takes the inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::EvaluateImpl(), which populates WorkPiece::outputs using the input arguments (passed to WorkPiece::EvaluateImpl()).  This function then checks WorkPiece::outputs, which much match WorkPiece::numOutputs and WorkPiece::outputTypes if they are specified.
	 
	 This function builds a reference vector to the inputs and calls WorkPiece::Evaluate(ref_vector<boost::any> const& ins)
	 @param[in] ins A vector of inputs 
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      std::vector<boost::any> Evaluate(std::vector<boost::any> const& ins);
      
      /// Evaluate this muq::Modeling::WorkPiece using references to the inputs
      /**
	 This function takes the references to the inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::EvaluateImpl(), which populates WorkPiece::outputs using the input arguments (passed to WorkPiece::EvaluateImpl()).  This function then checks WorkPiece::outputs, which much match WorkPiece::numOutputs and WorkPiece::outputTypes if they are specified.
	 
	 References are used for efficiency in the muq::Modeling::WorkGraph.
	 @param[in] ins A vector of references to the inputs 
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      std::vector<boost::any> Evaluate(ref_vector<boost::any> const& ins);
      
      /// Evaluate this muq::Modeling::WorkPiece in the case that there are no inputs
      /**
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      std::vector<boost::any> Evaluate();
      
      /// Evalaute this muq::Modeling::WorkPiece using multiple arguments
      /**
	 This function allows the user to call WorkPiece::Evaluate without first creating a vector of inputs.  Instead, the user calls WorkPiece::Evaluate with multiple arguments (if specified, the number of arguments must match the number of inputs) and this function creates the input vector.
	 @param[in] args The inputs (may be more than one)
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      template<typename... Args>			
	std::vector<boost::any> Evaluate(Args... args) {
	
	// we have new outputs
	outputs.clear();
	
	// create the reference input vector
	ref_vector<boost::any> inputs;
	inputs.reserve(numInputs<0? 0 : numInputs);
	
	// begin calling the EvaluateMulti with the first input
	return EvaluateRecursive(inputs, args...);
      }

      /// Get the (unique) name of this work piece
      std::string Name() const;
      
      /// The number of inputs
      const int numInputs;
	
      /// The number of outputs
      const int numOutputs;
      
    protected:
      
      /// The outputs
      /**
	 The outputs of this muq::Modeling::WorkPiece are filled by WorkPiece::EvaluateImpl().  If the number of outputs is specified (i.e., WorkPiece::numOutputs is not -1) then WorkPiece::Evaluate() checks to make sure the size of this vector is equal to WorkPiece::numOutputs after calling WorkPiece::EvaluateImpl().  If the output types are specified (i.e., WorkPiece::outputTypes is not an empty vector) then WorkPiece::Evaluate() checks that the output types match WorkPiece::outputTypes after calling WorkPiece::EvaluateImpl().
      */
      std::vector<boost::any> outputs = std::vector<boost::any>(0);
	
    private:
      
      /// The input types
      /**
	 Each element specifies the type of the corresponding input.  This vector must have the same number of elements as WorkPiece::numInputs or it is empty (default), which indicates that the input types are variable.
      */
      std::vector<std::string> inputTypes = std::vector<std::string>(0);
      
      /// The output types
      /**
	 Each element specifies the type of the corresponding output.  This vector must have the same number of elements as WorkPiece::numOutputs or it is empty (default), which indicates that the output types are variable.
      */
      std::vector<std::string> outputTypes = std::vector<std::string>(0);
      
      /// User-implemented function that determines the behavior of this muq::Modeling::WorkPiece
      /**
	 This function determines how the WorkPiece::inputs determine WorkPiece::outputs.  Must be implemented by a child.
	 
	 WorkPiece::Evaluate() calls this function after checking the inputs and storing them in WorkPiece::inputs.  This function populates WorkPiece::outputs, the outputs of this muq::Modeling::WorkPiece.  WorkPiece::Evaluate() checks the outputs after calling this function.
      */
      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) = 0;
      
      /// Creates WorkPiece::inputs when the WorkPiece::Evaluate is called with multiple arguments
      /**
	 @param[in] inputNum The current input number (the \f$i^{th}\f$ input)
	 @param[in] in The input corresponding to the \f$i^{th}\f$ input
	 @param[in] args The remaining (greater than \f$i\f$) inputs
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      template<typename ith, typename... Args>		       
	std::vector<boost::any> EvaluateRecursive(ref_vector<boost::any> &inputs, ith const& in, Args... args) {
	
	const int inputNum = inputs.size();
	
	// we have not yet put all of the inputs into the map, the ith should be less than the total number
	assert(numInputs<0 || inputNum+1<numInputs);
	
	if( inputTypes.size()>0 ) { // if we know the input types
	  // make sure the index is valid
	  assert(inputNum+1<inputTypes.size());
	  
	  // make sure the type is correct
	  assert(inputTypes[inputNum].compare(typeid(in).name())==0);
	}
	
	// add the last input to the input vector
	const boost::any in_any(in);
	inputs.push_back(std::cref(in_any));
	
	// call with EvaluateMulti with the remaining inputs
	return EvaluateRecursive(inputs, args...);
      }								
      
      /// Creates WorkPiece::inputs when the WorkPiece::Evaluate is called with multiple arguments
      /**
	 @param[in] inputNum The current input number (the last input, should match WorkPiece::numInputs-1 if it is specfied)
	 @param[in] in The input corresponding to the last input
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      template<typename last>			
	std::vector<boost::any> EvaluateRecursive(ref_vector<boost::any> &inputs, last const& in) {
	
	const int inputNum = inputs.size();
	
	// this is the last input, the last one should equal the total number of inputs
	assert(numInputs<0 || inputNum+1==numInputs);
	
	if( inputTypes.size()>0 ) { // if we know the input types
	  // make sure the index is valid
	  assert(inputNum+1==inputTypes.size());
	  
	  // make sure the type is correct
	  assert(inputTypes[inputNum].compare(typeid(in).name())==0);
	}
	
	// add the last input to the input vector
	const boost::any in_any(in);
	inputs.push_back(std::cref(in_any));
	
	return Evaluate(inputs);
      }

      /// Set the ID number, must be called by the constructor
      unsigned int SetID();

      /// A unique ID number assigned by the constructor
      const unsigned int id;
    }; // class WorkPiece
  } // namespace Modeling
} // namespace muq

#endif
