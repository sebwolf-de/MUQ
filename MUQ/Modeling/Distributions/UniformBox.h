#ifndef UNIFORMBOX_H_
#define UNIFORMBOX_H_

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"
#include "MUQ/Utilities/VariadicMacros.h"

#include "MUQ/Modeling/Distributions/Distribution.h"


namespace muq {
  namespace Modeling {

      /** Defines a normalized uniform distribution over a bounded rectangular domain. */
    class UniformBox : public Distribution {
    public:

      virtual ~UniformBox() = default;

      /// Create a \f$d\f$-dimensional uniform distribution
      /**
	 @param[in] bounds A \f$d\f$-dimensional vector, each entry is a pair---first: lower bound in the dimension, second: upper bound in that dimension
       */
      UniformBox(std::vector<std::pair<double, double> > const& bounds);

      /// Create a \f$d\f$-dimensional uniform distribution
      /**
	 @param[in] args The upper/lower bound pairs (may be more than one)
       */
      template<typename... Args>
      inline UniformBox(Args... args) : UniformBox(CreateBounds(args...)) {}
      
    private:

      /// Evaluate the log-density
      /**
	 Inputs:
	 <ol>
	 <li> The state \f$x\f$
	 </ol>
	 \return The log-density (either 1 or negative infinity)
       */
      virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) override;

      /// Sample the distribution
      virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs) override;

      STATIC_VARIADIC_TO_VECTOR(CreateBounds, (std::pair<double, double>), (std::vector<std::pair<double,double>>))
      inline static  std::vector<std::pair<double,double>> CreateBounds(std::vector<std::pair<double,double>>& bounds){return bounds;};
         
      /// Compute the bounds on the uniform box
      /**
	 @param bounds A vector of the bounds for each dimension---first: lower bound, second: upper bound
	 \return A vector of the first: lower bound, second: length for the box
       */
      static std::vector<std::pair<double, double> > ComputeBounds(std::vector<std::pair<double, double> > const& bounds);

      /// A vector that discribes the bounding box.
      /**
	 A vector that discribes the bounding box.  For each entry \f$i\f$, first: lower bound, second: length of the box for dimension \f$i\f$.  For example, \f$x \in [-2, 4]\f$ would be a pair with first entry \f$-2\f$ and second entry \f$6\f$.
       */
      const std::vector<std::pair<double, double> > bounds;

      /// The muq::Utilities::AnyAlgebra
      std::shared_ptr<muq::Utilities::AnyAlgebra> algebra;
    };
  } // namespace Modeling
} // namespace muq

#endif
