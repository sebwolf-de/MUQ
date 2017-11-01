
#ifndef REGRESSION_H_
#define REGRESSION_H_

#include <Eigen/QR>

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/AnyAlgebra.h"

#include "MUQ/Approximation/Regression/MultiIndex.h"
#include "MUQ/Approximation/Regression/Polynomial.h"

namespace muq {
  namespace Approximation {
    class Regression : public muq::Modeling::WorkPiece {
    public:

      /// Which polynomial basis should we use?
      enum PolynomialBasis {
	/// A monomial basis
	/**
	   May be poorly conditioned
	 */
	MonomialBasis,

	/// A Hermite polynomial
	/**
	   Orthogonal polynomial that may be unstable at high orders
	 */
	HermiteBasis,

	/// A Legendre polynomial
	/**
	   Orthogoal polynomial
	 */
	LegendreBasis
      };

      /**
	 @param[in] order The order of the polynomial regression
	 @param[in] basis The type of polynomial basis to use (defaults to Legendre)
       */
      Regression(unsigned int const order, Regression::PolynomialBasis const& basis = Regression::PolynomialBasis::LegendreBasis);

      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
	 @param[in] center The center of the inputs (used to recenter the inputs)
       */
      template<typename data>
	inline void Fit(std::vector<data> xs, std::vector<data> const& ys, boost::any const& center) {
	assert(xs.size()>0);
	assert(xs.size()==ys.size());

	// set the current center
	currentCenter = center;

	// center the input points
	CenterPoints(xs);
	
	// Compute basis coefficients
	ComputeCoefficients(xs, ys);
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

      int NumInterpolationPoints() const;
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// Compute the coefficients for the basis functions
      /**
	 Given points, data, and a center compute the coefficients on the basis (inner product with the basis evalautes the local polynomial).  The data can be a have multiple outputs (e.g., fitting more than one polynomial at once), which leads to more than one set of basis coefficents.
	 @param[in] xs The points 
	 @param[in] ys The output at each point 
      */
      template<typename data>
	inline void ComputeCoefficients(std::vector<data> const& xs, std::vector<data> const& ys) {
	assert(xs.size()==ys.size());

	// initalize the multi-index
	multi = std::make_shared<MultiIndex>(algebra->VectorDimensionBase(xs[0]), order);

	// check to make sure we have more than the number of points required to interpolate
	const unsigned int interp = NumInterpolationPoints();
	if( xs.size()<interp ) {
	  std::cerr << std::endl << "ERROR: Regression requires " << interp << " points to interpolate but only " << xs.size() << " are given." << std::endl;
	  std::cerr << "\tTry fitting the regression with at least " << interp+1 << " points." << std::endl << std::endl;
	  assert(xs.size()>NumInterpolationPoints());
	}

	// create the Vandermonde matrix and the rhs
	Eigen::MatrixXd vand = VandermondeMatrix(xs);
	const Eigen::MatrixXd rhs = ComputeCoefficientsRHS(vand, ys);
	vand = vand.transpose()*vand;

	// make the solver to do the regression
	auto solver = vand.colPivHouseholderQr();

	// comptue the coefficients
	coeff = solver.solve(rhs).transpose();
      }

      /// Compute the right hand side given data to compute the polynomial coefficients
      /**
	 @param[in] vand The Vandermonde matrix
      	 @param[in] ys_data The output at each point 
      */
      template<typename data>
	inline Eigen::MatrixXd ComputeCoefficientsRHS(Eigen::MatrixXd const& vand, std::vector<data> const& ys_data) const {
	// the dimension
	assert(ys_data.size()>0);
	const unsigned int dim = algebra->VectorDimensionBase(ys_data[0]);

	// initialize space for the data
	Eigen::MatrixXd ys = Eigen::MatrixXd::Constant(ys_data.size(), dim, std::numeric_limits<double>::quiet_NaN());

	// copy the data into an Eigen type
	for( unsigned int i=0; i<ys_data.size(); ++i ) {
	  for( unsigned int j=0; j<dim; ++j ) {
	    ys(i, j) = boost::any_cast<double const>(algebra->AccessElementBase(j, ys_data[i])); 
	  }
	}

	// apply the Vandermonde matrix
	return vand.transpose()*ys;
      }

      /// Create the Vandermonde matrix
      /**
	 @param[in] xs The points 
      */
      template<typename data>
	Eigen::MatrixXd VandermondeMatrix(std::vector<data> const& xs) const {
	assert(multi);
	
	// the number of points and the number of terms
	const unsigned int N = xs.size();
	const unsigned int M = multi->Size();
	assert(N>0);

	// initialize the matrix
	Eigen::MatrixXd vand = Eigen::MatrixXd::Ones(N, M);

	// each term is built by evaluating the polynomial basis
	for( unsigned int pt=0; pt<N; ++pt ) { // loop through the points
	  for( unsigned int i=0; i<M; ++i ) { // loop through the terms
	    // get the multi-index
	    const Eigen::VectorXi alpha = boost::any_cast<Eigen::VectorXi const&>(multi->Evaluate(i) [0]);

	    // get the point
	    const data& pnt = xs[pt];
	    assert(alpha.size()==algebra->VectorDimensionBase(pnt));

	    // each term is a product of 1D variables
	    for( unsigned int v=0; v<alpha.size(); ++v ) {
	      // the point where we are evaluating the polynomial
	      const double x = boost::any_cast<double const>(algebra->AccessElementBase(v, pnt));

	      // evaluate the polynomial
	      vand(pt, i) *= boost::any_cast<double const>(poly->Evaluate((unsigned int)alpha(v), x) [0]);
	    }
	  }
	}

	return vand;
      }
      
      /// Center the input points
      template<typename data>
	inline void CenterPoints(std::vector<data>& xs) {
	// reset the current radius
	currentRadius = 0.0;

	// is the center zero?
	const bool zeroCenter = algebra->IsZeroBase(currentCenter);
	
	// loop through all of the input points
	for( auto it=xs.begin(); it!=xs.end(); ++it ) {
	  if( !zeroCenter ) {
	    // recenter the the point
	    const boost::any vec = *it;
	    *it = boost::any_cast<data const&>(algebra->SubtractBase(std::reference_wrapper<boost::any const>(vec), std::reference_wrapper<boost::any const>(currentCenter)));
	  }
	  
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

      /// The order of the regression
      const unsigned int order;

      /// An muq::Modeling::AnyAlgebra to do the algebric manipulations
      std::shared_ptr<muq::Modeling::AnyAlgebra> algebra;

      /// The multi-index to so we know the order of each term
      std::shared_ptr<MultiIndex> multi;

      /// The polynomial basis (in one variable) used to compute the Vandermonde matrix
      std::shared_ptr<Polynomial> poly;

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

      /// Coeffients for the polynomial basis
      Eigen::MatrixXd coeff;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
