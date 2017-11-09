#ifndef UNIFORM_H_
#define UNIFORM_H_

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace Modeling {
    class Uniform : public Distribution {
    public:

      virtual ~Uniform();

      /// Create a \f$d\f$-dimensional uniform distribution
      /**
	 @param[in] bounds A \f$d\f$-dimensional vector, each entry is a pair---first: lower bound in the dimension, second: upper bound in that dimension
       */
      Uniform(std::vector<std::pair<double, double> > const& bounds);

      /// Create a \f$d\f$-dimensional uniform distribution
      /**
	 @param[in] args The upper/lower bound pairs (may be more than one)
       */
      template<typename... Args>
	inline Uniform(Args... args) : Uniform(CreateBounds(args...)) {}
      
    private:

      /// Evaluate the log-density
      /**
	 Inputs:
	 <ol>
	 <li> The state \f$x\f$
	 </ol>
	 \return The log-density (either 1 or negative infinity)
       */
      virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) const override;

      /// Sample the distribution
      virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs) override;

      /// Create a vector of the upper/lower bound pairs
      /**
	 @param[in] args The upper/lower bound pairs (may be more than one)
       */
      template<typename... Args>
	inline static std::vector<std::pair<double, double> > CreateBounds(Args... args) {
	// create an empty vector
	std::vector<std::pair<double, double> > bounds(0);

	// return the upper/lower bound pairs
	return CreateBounds(bounds, args...);
      }

      /// Create a vector of the upper/lower bound pairs
      /**
	 @param[out] bounds The uppwer/lower bound pairs so far
	 @param[in] ith The \f$i^{th}\f$ upper/lower boundd pair
	 @param[in] args The remaining upper/lower bound pairs (may be more than one)
       */
      template<typename... Args>
	inline static std::vector<std::pair<double, double> > CreateBounds(std::vector<std::pair<double, double> >& bounds, std::pair<double, double> ith, Args... args) {
	// add the ith pair
	bounds.push_back(ith);

	// add the remaining pairs
	return CreateBounds(bounds, args...);
      }

      /// Create a vector of the upper/lower bound pairs
      /**
	 @param[out] bounds The uppwer/lower bound pairs so far
	 @param[in] last The last upper/lower boundd pair
       */
      static std::vector<std::pair<double, double> > CreateBounds(std::vector<std::pair<double, double> >& bounds, std::pair<double, double> last);

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
