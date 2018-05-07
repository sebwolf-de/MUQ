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


      /// Are we specifying the mean, covariance matrix, or precision matrix
      enum ExpectedInputs {
        LocationOnly, /// Only the location can be specified.  All other inputs are fixed apriori.
        LocationMean, /// Only the location and mean need to be specified
        All,          /// Every evaluation requires specifying all inputs
        Flexible      /// Input type may change with every evaluation
      };

      /// Construct a Gaussian with scaled identity covariance/precision
      /**
	 @param[in] obj Either the mean, covariance, or the precision (depending on the second paameter)
	 @param[in] mode Are we specifying mean, covariance, or precision (defaults to mean)
       */
      Gaussian(boost::any               const& obj = 0.0,
               Gaussian::Mode           const  mode = Gaussian::Mode::Mean,
               Gaussian::ExpectedInputs const  expected = Gaussian::ExpectedInputs::LocationOnly);

      /// Construct a Gaussian with scaled identity covariance/precision
      /**
	 @param[in] mean The mean
	 @param[in] obj Either the covariance or the precision (depending on the second paameter)
	 @param[in] mode Are we specifying mean, covariance, or precision (defaults to covariance)
       */
      Gaussian(boost::any               const& mean,
               boost::any               const& obj,
               Gaussian::Mode           const  mode = Gaussian::Mode::Covariance,
               Gaussian::ExpectedInputs const  expected = Gaussian::ExpectedInputs::LocationOnly);

      ~Gaussian();

      /// Get the dimension of this Gaussian
      /**
	 \return The dimension
       */
      unsigned int Dimension() const;

      /// Get the covariance
      /**
	 \return The covariance
       */
      boost::any GetCovariance() const;

      /// Set the covariance matrix
      /**
	 @param[in] newcov The new covariance
       */
      void SetCovariance(boost::any const& newcov);

      /// Set the precision matrix
      /**
	 @param[in] newprec The new precision
       */
      void SetPrecision(boost::any const& newprec);

      void SetInputTypes(Gaussian::ExpectedInputs const  expected);
      
    private:

      /// Compute the distribution's scaling constant
      /**
	 \return Scaling constant
       */
      void ComputeScalingConstant();

      static boost::any SaveCovPrec(boost::any const& in);

      /// Reset the hyperparameters
      /**
	 @param[in] hyperparas A list of hyperparameters mean and/or covariance/precision
       */
      void ResetHyperparameters(ref_vector<boost::any> const& hyperparas);

      /// Implement the log-density for a Gaussian distribution
      /**
	 Inputs:
	 <ol>
	 <li> The state \f$x\f$
	 </ol>
	 \return The log-density
       */
      virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) override;

      /// Sample the distribution
      virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs) override;

      /// The muq::Utilities::AnyAlgebra
      std::shared_ptr<muq::Utilities::AnyAlgebra> algebra = std::make_shared<muq::Utilities::AnyAlgebra>();

      /// Have we specified the covariance or the precision
      Gaussian::Mode mode;

      /// The dimension
      const unsigned int dim;

      /// The mean of the distribution
      boost::optional<boost::any> mean;

      /// The covariance
      boost::optional<boost::any> cov;

      /// The square root of the covariance
      boost::optional<boost::any> covSqrt;

      /// The precision
      boost::optional<boost::any> prec;

      /// The square root of the precision
      boost::optional<boost::any> precSqrt;

      /// The scaling constant for the density
      double scalingConstant;

    };
  } // namespace Modeling
} // namespace muq

#endif
