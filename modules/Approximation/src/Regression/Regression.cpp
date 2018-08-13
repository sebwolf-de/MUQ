#include "MUQ/Approximation/Regression/Regression.h"

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;

Regression::Regression(unsigned int const order, std::string const& polyName) : WorkPiece(), order(order) {
  poly = IndexedScalarBasis::Construct(polyName);
}

void Regression::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // if there are no points ... just return with an empty outputs
  if(inputs.size()==0) { return; }

  // make sure we can compute the Vandermonde matrix
  assert(multi);

  std::vector<Eigen::VectorXd> centered(inputs.size());
  assert(currentRadius>0.0);
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    // center and normalize
    centered[i] = (boost::any_cast<Eigen::VectorXd const&>(inputs[i])-currentCenter)/currentRadius;
  }

  // get the Vandermonde matrix of the inputs
  const Eigen::MatrixXd vand = VandermondeMatrix(centered);
  assert(coeff.cols()==vand.cols());

  // compute the regression polynomial
  outputs.resize(1);
  outputs[0] = (Eigen::MatrixXd)(coeff*vand.transpose());
}

int Regression::NumInterpolationPoints() const {
  if( multi ) {
    return multi->Size();
  }

  std::cerr << std::endl << std::endl << "ERROR: Not able to compute the number of points required for interpolation" <<
    std::endl << "\tPolynomialRegressor.cpp NumInterpolationPoints()" << std::endl;
  assert(false);
  
  return -1;
}

void Regression::Fit(std::vector<Eigen::VectorXd> xs, std::vector<Eigen::VectorXd> const& ys, Eigen::VectorXd const& center) {
  assert(xs.size()>0);
  assert(xs.size()==ys.size());
  
  // set the current center
  currentCenter = center;
  
  // center the input points
  CenterPoints(xs);
  
  // Compute basis coefficients
  ComputeCoefficients(xs, ys);
}

void Regression::Fit(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys) {
  assert(xs.size()>0);
  assert(xs.size()==ys.size());
  
  // preform the fit with zero center
  Fit(xs, ys, Eigen::VectorXd::Zero(xs[0].size()));
}

void Regression::ComputeCoefficients(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys) {
  assert(xs.size()==ys.size());
  
  // initalize the multi-index
  multi = MultiIndexFactory::CreateTotalOrder(xs[0].size(), order);
  
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

Eigen::MatrixXd Regression::ComputeCoefficientsRHS(Eigen::MatrixXd const& vand, std::vector<Eigen::VectorXd> const& ys_data) const {
  // the dimension
  assert(ys_data.size()>0);
  const unsigned int dim = ys_data[0].size();
  
  // initialize space for the data
  Eigen::MatrixXd ys = Eigen::MatrixXd::Constant(ys_data.size(), dim, std::numeric_limits<double>::quiet_NaN());
  
  // copy the data into an Eigen type
  for( unsigned int i=0; i<ys_data.size(); ++i ) {
    ys.row(i) = ys_data[i];
  }
  
  // apply the Vandermonde matrix
  return vand.transpose()*ys;
}

Eigen::MatrixXd Regression::VandermondeMatrix(std::vector<Eigen::VectorXd> const& xs) const {
  assert(multi);
  
  // the number of points and the number of terms
  const unsigned int N = xs.size();
  const unsigned int M = multi->Size();
  assert(N>0);
  
  // initialize the matrix
  Eigen::MatrixXd vand = Eigen::MatrixXd::Ones(N, M);
  
  // each term is built by evaluating the polynomial basis
  for( unsigned int i=0; i<M; ++i ) { // loop through the terms
    // get the multi-index
    const Eigen::RowVectorXi& alpha = multi->at(i)->GetVector();
    
    for( unsigned int pt=0; pt<N; ++pt ) { // loop through the points
      // get the point
      const Eigen::VectorXd& pnt = xs[pt];
      assert(alpha.size()==pnt.size());
      
      // each term is a product of 1D variables
      for( unsigned int v=0; v<alpha.size(); ++v ) {
	// evaluate the polynomial
	vand(pt, i) *= boost::any_cast<double const>(poly->Evaluate((unsigned int)alpha(v), pnt[v]) [0]);
      }
    }
  }
  
  return vand;
}

void Regression::CenterPoints(std::vector<Eigen::VectorXd>& xs) {
  // reset the current radius
  currentRadius = 0.0;
  
  // is the center zero?
  const bool zeroCenter = currentCenter.norm()<std::numeric_limits<double>::epsilon();
  
  // loop through all of the input points
  for( auto it=xs.begin(); it!=xs.end(); ++it ) {
    if( !zeroCenter ) {
      // recenter the the point
      *it -= currentCenter;
    }
    
    // set the radius to the largest distance from the center
    currentRadius = std::max(currentRadius, it->norm());
  }
  
  // loop through all of the input points to normalize by the radius
  for( auto it=xs.begin(); it!=xs.end(); ++it ) {
    *it /= currentRadius;
  }
}
