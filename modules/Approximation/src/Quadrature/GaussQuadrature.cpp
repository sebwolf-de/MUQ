#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"

using namespace muq::Approximation;

GaussQuadrature::GaussQuadrature() {}

// Will this constructor do anything else?
// Is this how we can define polyOrder?
GaussQuadrature::GaussQuadrature(std::shared_ptr<Polynomial> poly_in)
                                : poly(poly_in), polyOrder(poly_in->order()) {}

void GaussQuadrature::Calculate() {

  // Create diagonal and subdiagonal vectors
  Eigen::VectorXd diag = Eigen::VectorXd::Zero(polyOrder);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Zero(polyOrder);

  for (unsigned int i=1; i<polyOrder+1; i++) {

    // Need to connect to the a,b, and c functions for the specic polynomial
    double alpha_i = -b(i)/a(i);
    double beta_i = std::sqrt(c(i+1)/(a(i)*a(i+1)));

    // Diagonal entry of J
    diag(i-1) = alpha_i;

    // Off diagonal entries of J
    if (i < polyOrder)
      J(i-1) = beta_i;

  }

  Eigen::SelfAdjointEigenSolver<MatrixXd> es;
  es.computeFromTridiagonal(diag, subdiag);

  // Set gauss points
  // Will the eigenvalues function return an Eigen::VectorXd?
  // Can you assign a VectorXd this way?
  gauss_pts = es.eigenvalues();

  // Get mu_0 value (integral of weighting function)
  // Where do we get this from?
  double mu_0 = poly->mu;

  for (unsigned int i=0; i<polyOrder; i++) {

    // Check if these eigenvectors have magnitude 1.0
    Eigen::VectorXd eigen_vector_i = es.eigenvectors().col(i);

    // Set gauss weights
    gauss_wts(i) = mu_0*eigen_vector_i(0)**2;

  }

}

Eigen::VectorXd const& GaussQuadrature::points() const {

  return gauss_pts;

}

Eigen::VectorXd const& GaussQuadrature::weights() const {

  return gauss_wts;

}
