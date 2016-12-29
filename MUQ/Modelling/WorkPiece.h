#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<cassert>
#include<memory>

#include "boost/any.hpp"

namespace muq {
  namespace Modelling {
    /// Base class for MUQ's modeling envronment
    class WorkPiece {
    public:

      WorkPiece();

      WorkPiece(unsigned int const num, bool const fixInput = true);
            
      WorkPiece(unsigned int const numIns, unsigned int const numOuts);

      WorkPiece(std::vector<std::string> const& types, bool const fixInput = true);

      WorkPiece(std::vector<std::string> const& types, unsigned int const num, bool const fixInputType = true);

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
