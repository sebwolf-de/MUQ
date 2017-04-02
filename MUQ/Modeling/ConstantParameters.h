#ifndef CONSTANTPARAMETERS_H_
#define CONSTANTPARAMETERS_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    /// A muq::Modeling::WorkPiece with no inputs and known, constant outputs
    class ConstantParameters : public WorkPiece {
    public:

      /// Create a muq::Modeling::ConstantParameters with the outputs given in a vector
      /**
	 @param[in] outs The outputs
       */
      ConstantParameters(std::vector<boost::any> const& outs);

      /// Create a muq::Modeling::ConstantParameters with the outputs given indivually
      /**
	 @param[in] in The \f$i^{th}]f$ input
	 @param[in] args The \f$[i+1, N\f$] outputs
       */
      template<typename ith, typename... Args>
	ConstantParameters(ith const& in, Args... args) : ConstantParameters(args...) {
	// add this output to the begining of the output vector
	outputs.insert(outputs.begin(), in);

	// set the number of outputs
	numOutputs = outputs.size();

	// set the output types
	outputTypes = Types(Types(outputs));
      }

      /// Create a muq::Modeling::ConstantParameters with the outputs given indivually
      /**
	 @param[in] last the last output
      */
      template<typename last>
	ConstantParameters(last const& in) : WorkPiece(0, WorkPiece::Fix::Inputs) {	
	// the outputs will never change so we should not clear them
	clearOutputs = false;

	// add this output to the begining of the output vector
	outputs.insert(outputs.begin(), in);

	// set the number of outputs
	numOutputs = outputs.size();

	// set the output types
	outputTypes = Types(Types(outputs));
      }

    private:
      /// The outputs are already set and not cleared so don't do anything
      /**
	 @param[in] inputs An empty vector of inputs
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// The outputs, calling Evaluate will always return these values
      std::vector<boost::any> outs;

    };
  } // namespace Modeling
} // namespace muq

#endif
