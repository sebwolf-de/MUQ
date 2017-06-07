#ifndef ANYALGEBRA_H_
#define ANYALGEBRA_H_

#include <iostream>
#include <assert.h>

#include "boost/none.hpp"
#include "boost/any.hpp"

#include <Eigen/Core>

namespace muq {
  namespace Modeling {
    /// Implement a generic way to do algebric operations on boost::any's
    class AnyAlgebra {
    public:

      /// Default constructor
      AnyAlgebra();

      /// Compute an identity object for boost::any
      /**
	 For example, if the underlying type is a double this would return 1.0, if the underlying type is an Eigen::VectorXd this would return Eigen::MatrixXd::Zero(N).
	 @param[in] in The input --- its type and size determines the return type and size
	 \return An identity of some type
       */
      boost::any IdentityBase(std::reference_wrapper<const boost::any> const& in) const;
      
      /// Add two objects together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      boost::any AddBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const;

      /// Multiply two objects 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1) --- order matters!
       */
      boost::any MultiplyBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const;

    private:

      /// Compute an identity object for boost::any
      /**
	 MUQ automatically checks for some common input types.  However, the user may need to overload this function for special types.
	 @param[in] in The input --- its type and size determines the return type and size
	 \return An identity of some type
       */
      virtual boost::any Identity(std::reference_wrapper<const boost::any> const& in) const;

      /// Add two objects together
      /**
	 MUQ automatically checks for some common pairs.  However, the user may need to overload this function for special types.
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      virtual boost::any Add(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const;

      /// Multiply two objects 
      /**
	 MUQ automatically checks for some common pairs.  However, the user may need to overload this function for special types.
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1) --- order matters!
       */
      boost::any Multiply(std::reference_wrapper<const boost::any> const& in, std::reference_wrapper<const boost::any> const& out) const;

      const std::string eigenVecType = typeid(Eigen::VectorXd).name();

      const std::string eigenMatType = typeid(Eigen::MatrixXd).name();
      
      const std::string doubleType = typeid(double).name();

      const std::string noneType = typeid(boost::none).name();
    };
  } // namespace Modeling
} // namespace muq

#endif
