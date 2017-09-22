#ifndef REGRESSION_H_
#define REGRESSION_H_

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/AnyAlgebra.h"

namespace muq {
  namespace Approximation {
    class Regression : public muq::Modeling::WorkPiece {
    public:

      Regression();

      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
	 @param[in] center The center of the inputs (used to recenter the inputs)
       */
      template<typename data>
	inline void Fit(std::vector<data> xs, std::vector<data> const& ys, data const& center) {
	assert(xs.size()>0);
	assert(xs.size()==ys.size());

	// set the current center
	currentCenter = center;

	// center the input points
	CenterPoints(xs);

	// Compute basis coefficients
	ComputeCoefficients(xs, ys);
	
	std::cout << "INSIDE FIT" << std::endl;
      }
      
      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
       */
      template<typename data>
	inline void Fit(std::vector<data> const& xs, std::vector<data> const& ys) {
	assert(xs.size()>0);
	assert(xs.size()==ys.size());
	
	// get the zero of this vector type
	const unsigned int size = algebra->VectorDimensionBase(xs[0]);
	const boost::any& zero = algebra->ZeroVectorBase(typeid(data).name(), size);

	// preform the fit with zero center
	Fit<data>(xs, ys, boost::any_cast<data const&>(zero));
      }
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// Compute the coefficients for the basis functions
      /**
	 Given points, data, and a center compute the coefficients on the basis (inner product with the basis evalautes the local polynomial).  The data can be a have multiple outputs (e.g., fitting more than one polynomial at once), which leads to more than one set of basis coefficents.
	 @param[in] xs The points 
	 @param[in] ys The output at each point 
      */
      template<typename data>
	inline void ComputeCoefficients(std::vector<data> const& xs, std::vector<data> const& ys) const {
	assert(xs.size()==ys.size());

	// set the weights equal to one
	const Eigen::VectorXd weights = Eigen::VectorXd::Ones(ys.size());

	// create the Vandermonde matrix
	const Eigen::MatrixXd vand = VandermondeMatrix(xs);

	std::cout << "HERE" << std::endl;
      }

      // Create the Vandermonde matrix
      /**
	 @param[in] xs The points 
      */
      template<typename data>
	Eigen::MatrixXd VandermondeMatrix(std::vector<data> const& xs) const {
	// the number of points
	const unsigned int N = xs.size();
	
	return Eigen::MatrixXd();
      }
      
      /// Center the input points
      template<typename data>
	inline void CenterPoints(std::vector<data>& xs) {
	// if the current center is zero, do nothing
	if( algebra->IsZeroBase(currentCenter) ) {
	  return;
	}

	// reset the current radius
	currentRadius = 0.0;

	// loop through all of the input points
	for( auto it=xs.begin(); it!=xs.end(); ++it ) {
	  const boost::any vec = *it;

	  // recenter the the point
	  *it = boost::any_cast<data const&>(algebra->SubtractBase(std::reference_wrapper<boost::any const>(vec), std::reference_wrapper<boost::any const>(currentCenter)));

	  // set the radius to the largest distance from the center
	  currentRadius = std::max(currentRadius, algebra->NormBase(*it));
	}

	// loop through all of the input points to normalize by the radius
	const boost::any normalize = 1.0/currentRadius;
	for( auto it=xs.begin(); it!=xs.end(); ++it ) {
	  const boost::any vec = *it;
	  *it = boost::any_cast<data const&>(algebra->MultiplyBase(normalize, vec));
	}
      }

      /// An muq::Modeling::AnyAlgebra to do the algebric manipulations
      std::shared_ptr<muq::Modeling::AnyAlgebra> algebra;

      /// Current center of the inputs
      /**
	 Defaults to boost::none.
       */
      boost::any currentCenter = boost::none;

      /// Current radius of inputs
      /**
	 Defaults to zero.
       */
      double currentRadius = 0.0;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
