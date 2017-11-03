#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/AnyAlgebra.h"

namespace muq {
  namespace Modeling {
    class Gaussian : public Distribution {
    public:

      /// Are we specifying the covariance or the precision matrix
      enum Mode {
	/// We are specifying the covariance
	Covariance,

	/// We are specifying the precision
	Precision
      };

      /// Construct a zero mean Gaussian with scaled identity covariance/precision
      /**
	 @param[in] dim The dimension of the state (defaults to 1)
	 @param[in] cov_prec The covariance or the precision; depending on the mode
	 @param[in] mode Are we specifying a scaled identity covariance or precision (defaults to covariance)
       */
      Gaussian(unsigned int const dim = 1, double const cov_prec = 1.0, Gaussian::Mode const mode = Gaussian::Mode::Covariance);

      ~Gaussian();
      
    private:

      /// Implement the log-density for a Gaussian distribution
      /**
	 Inputs:
	 <ol>
	 <li> The state \f$x\f$
	 </ol>
	 \return The log-density 
       */
      virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) const override;

      /// Sample the distribution
      virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs) const override;

      /// The muq::Modeling::AnyAlgebra
      std::shared_ptr<AnyAlgebra> algebra;

      /// The dimension 
      const unsigned int dim;

      /// The covariance
      const double cov;
    };
  } // namespace Modeling
} // namespace muq

#endif
