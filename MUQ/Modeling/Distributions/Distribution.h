#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    class Distribution : public WorkPiece {
    public:

      /// Are we evaluting the log-density or sampling?
      enum Mode {
	/// Evaluate the log-density
	EvaluateLogDensity,

	/// Sample the distribution
	SampleDistribution
      };

      virtual ~Distribution();

      /// Create a muq::Modeling::Distribution with no fixed number of inputs and variable input/output types
      Distribution();

      /// Evaluate the log-density
      /**
	 If known, the log-density should be implemented by a child in the LogDensityImpl class.
	 @param[in] inputs the vector of inputs to the log-density
	 \return The log density
       */
      virtual double LogDensity(ref_vector<boost::any> const& inputs);

      /// Evaluate the log-density
      /**
	 This allows the user to evaluate the log-density without first creating a vector of inputs
	 @param[in] args The inputs (may be more than one)
	 \return The log density
       */
      template<typename... Args>
	inline double LogDensity(Args... args) {
	// create the reference input vector
	ref_vector<boost::any> inputs;
	inputs.reserve(numInputs<0? 0 : numInputs-1); // the first input is always whether we are evaluting the log-density or sampling

	// begin calling LogDensity recursively
	return LogDensity(inputs, args...);
      }

      /// Sample the distribution
      /**
	 Calls SampleImpl, the default behavior is to return boost::none
	 @param[in] inputs the vector of inputs to the log-density
	 \return A sample
       */
      boost::any Sample(ref_vector<boost::any> const& inputs);

      /// Sample the distribution with no inputs
      /**
	 Allows the user to call Sample without any inputs
	 \return A sample
       */
      boost::any Sample();

      /// Sample the distribution
      /**
	 This allows the sample the distribution without first creating a vector of inputs
	 @param[in] args The inputs (may be more than one)
	 \return A sample
       */
      template<typename... Args>
	inline boost::any Sample(Args... args) {
	// create the reference input vector
	ref_vector<boost::any> inputs;
	inputs.reserve(numInputs<0? 0 : numInputs-1); // the first input is always whether we are evaluting the log-density or sampling

	// begin calling Sample recursively
	return Sample(inputs, args...);
      }
      
    private:

      /// Implement the log-density
      /**
	 If known, the log-density should be implemented by a child.  If it is not overridden then the default behavior is to return negative infinity (-1.0*std::numeric_limits<double>::infinity()).
	 @param[in] inputs the vector of inputs to the log-density
	 \return The log density
       */
      virtual double LogDensityImpl(ref_vector<boost::any> const& inputs);

      /// Sample the distribution
      /**
	 Should be overwritten by a child.  The default behavior is to return boost::none
       */
      virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs);

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      /// Evaluate the log-density
      /**
	 This allows the user to evaluate the log-density without first creating a vector of inputs
	 @param[in] inputs The inputs so far
	 @param[in] in The \f$i^{th}\f$ input
	 @param[in] args The remaining inputs (may be more than one)
	 \return The log density
       */
      template<typename ith, typename... Args>
	inline double LogDensity(ref_vector<boost::any>& inputs, ith const& in, Args... args) {
	const int inputNum = inputs.size()+1; // add one, the first input is always whether we are evaluting the log-density or sampling
		
	// we have not yet put all of the inputs into the map, the ith should be less than the total number
	assert(numInputs<0 || inputNum<numInputs);

	// check the input type
	assert(CheckInputType(inputNum, typeid(in).name()));
	
	// add the last input to the input vector
	const boost::any in_any(in);
	inputs.push_back(std::cref(in_any));
	
	// call LogDensity recursively
	return LogDensity(inputs, args...);
      }

      /// Sample the distribution
      /**
	 This allows the user to sample the distribution without first creating a vector of inputs
	 @param[in] inputs The inputs so far
	 @param[in] in The \f$i^{th}\f$ input
	 @param[in] args The remaining inputs (may be more than one)
	 \return A sample
       */
      template<typename ith, typename... Args>
	inline boost::any Sample(ref_vector<boost::any>& inputs, ith const& in, Args... args) {
	const int inputNum = inputs.size()+1; // add one, the first input is always whether we are evaluting the log-density or sampling
		
	// we have not yet put all of the inputs into the map, the ith should be less than the total number
	assert(numInputs<0 || inputNum<numInputs);

	// check the input type
	assert(CheckInputType(inputNum, typeid(in).name()));
	
	// add the last input to the input vector
	const boost::any in_any(in);
	inputs.push_back(std::cref(in_any));
	
	// call Sample recursively
	return Sample(inputs, args...);
      }

      /// Add the mode to the begining of an input type vector
      /**
	 @param[in] types A vector of input types
	 \return The vector of input times with mode added to the begining
       */
      static std::vector<std::string> AddModeInput(std::vector<std::string> const& types);

      /// Add the mode to the begining of an input type map
      /**
	 @param[in] types A map from input number to input types
	 \return The map from input number to input types with the mode added as the first one
       */
      static std::map<unsigned int, std::string> AddModeInput(std::map<unsigned int, std::string> const& types);
    };
  } // namespace Modeling
} // namespace muq

#endif
