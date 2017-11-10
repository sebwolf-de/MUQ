#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace Modeling {
    class Gaussian : public Distribution {
    public:

      /// Are we specifying the mean, covariance matrix, or precision matrix
      enum Mode {
	/// We are specifying the mean
	Mean,

	/// We are specifying the covariance
	Covariance,

	/// We are specifying the precision
	Precision
      };

      /// Construct a Gaussian with scaled identity covariance/precision
      /**
	 @param[in] obj Either the mean, covariance, or the precision (depending on the second paameter)
	 @param[in] mode Are we specifying mean, covariance, or precision (defaults to mean)
       */
      Gaussian(boost::any const& obj = 0.0, Gaussian::Mode const mode = Gaussian::Mode::Mean);

      /// Construct a Gaussian with scaled identity covariance/precision
      /**
	 @param[in] mean The mean
	 @param[in] obj Either the covariance or the precision (depending on the second paameter)
	 @param[in] mode Are we specifying mean, covariance, or precision (defaults to covariance)
       */
      Gaussian(boost::any const& mean, boost::any const& obj, Gaussian::Mode const mode = Gaussian::Mode::Covariance);

      ~Gaussian();

      /// Get the dimension of this Gaussian
      /**
	 \return The dimension
       */
      unsigned int Dimension() const;
      
    private:

      static boost::any SaveCovPrec(boost::any const& in);

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
      virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs) override;

      /// The muq::Utilities::AnyAlgebra
      std::shared_ptr<muq::Utilities::AnyAlgebra> algebra = std::make_shared<muq::Utilities::AnyAlgebra>();

      /// Have we specified the covariance or the precision
      const Gaussian::Mode mode;

      /// The dimension 
      const unsigned int dim;

      /// The mean of the distribution
      const boost::optional<boost::any> mean;

      /// The covariance 
      const boost::optional<boost::any> cov;

      /// The square root of the covariance
      boost::optional<boost::any> covSqrt;

      /// The precision
      const boost::optional<boost::any> prec;

      /// The square root of the precision
      boost::optional<boost::any> precSqrt;

    };
  } // namespace Modeling
} // namespace muq

#endif
