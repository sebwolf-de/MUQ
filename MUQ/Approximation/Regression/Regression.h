#ifndef REGRESSION_H_
#define REGRESSION_H_

#include <Eigen/QR>

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"
#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

namespace muq {
  namespace Approximation {
    class Regression : public muq::Modeling::WorkPiece {
    public:
      
      /**
	 @param[in] order The order of the polynomial regression
	 @param[in] basis The type of polynomial basis to use (defaults to Legendre)
      */
      Regression(unsigned int const order, std::string const& basis = "Legendre");
      
      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
	 @param[in] center The center of the inputs (used to recenter the inputs)
      */
      void Fit(std::vector<Eigen::VectorXd> xs, std::vector<Eigen::VectorXd> const& ys, Eigen::VectorXd const& center);
      
      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
       */
      void Fit(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys);

      int NumInterpolationPoints() const;
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// Compute the coefficients for the basis functions
      /**
	 Given points, data, and a center compute the coefficients on the basis (inner product with the basis evalautes the local polynomial).  The data can be a have multiple outputs (e.g., fitting more than one polynomial at once), which leads to more than one set of basis coefficents.
	 @param[in] xs The points 
	 @param[in] ys The output at each point 
      */
      void ComputeCoefficients(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys); 

      /// Compute the right hand side given data to compute the polynomial coefficients
      /**
	 @param[in] vand The Vandermonde matrix
      	 @param[in] ys_data The output at each point 
      */
      Eigen::MatrixXd ComputeCoefficientsRHS(Eigen::MatrixXd const& vand, std::vector<Eigen::VectorXd> const& ys_data) const;

      /// Create the Vandermonde matrix
      /**
	 @param[in] xs The points 
	 \return The Vandermonde matrix
      */
      Eigen::MatrixXd VandermondeMatrix(std::vector<Eigen::VectorXd> const& xs) const;
      
      /// Center the input points
      void CenterPoints(std::vector<Eigen::VectorXd>& xs);

      /// The order of the regression
      const unsigned int order;

      /// The multi-index to so we know the order of each term
      std::shared_ptr<muq::Utilities::MultiIndexSet> multi;

      /// The polynomial basis (in one variable) used to compute the Vandermonde matrix
      std::shared_ptr<IndexedScalarBasis> poly;

      /// Current center of the inputs
      Eigen::VectorXd currentCenter;

      /// Current radius of inputs
      /**
	 Defaults to zero.
       */
      double currentRadius = 0.0;

      /// Coeffients for the polynomial basis
      Eigen::MatrixXd coeff;
    };
  } // namespace Approximation
} // namespace muq

#endif
