#include "MUQ/Approximation/Polynomials/GaussQuadrature.h"

using namespace muq::Approximation;

// void GaussQuadrature::Calculate(){
//
//   std::vector<Eigen::Triplet<double>> triple_list;
//
//   // Where do we get polyOrder?
//   // Loop to create triplet list
//   for (unsigned int i=0; i<polyOrder; i++){
//
//     // Need to connect to the a,b, and c functions for the specic polynomial
//     double alpha_i = -b(i)/c(i);
//     double beta_i = std::sqrt(c(i+1)/(a(i)*a(i+1)));
//
//     // Diagonal entry of J
//     triple_list.push_back(Eigen::Triplet<double>(i,i,alpha_i));
//
//     // Off diagonal entries of J
//     if (i < polyOrder-1){
//       triple_list.push_back(Eigen::Triplet<double>(i,i+1,beta_i));
//       triple_list.push_back(Eigen::Triplet<double>(i+1,i,beta_i));
//     }
//
//   }
//
//   // Create sparse matrix J
//   Eigen::SparseMatrix<double> J(polyOrder,polyOrder)
//
//   // Assign entries of the sparse matrix from triplet list
//   J.setFromTriplets(triple_list.begin(), triple_list.end());
//
//   // Set gauss points
//   gauss_pts =
//
//   // Set gauss weights
//   gauss_weights =
//
// }


GaussQuadrature::GaussQuadrature() {}

// Will this constructor do anything else?
GaussQuadrature::GaussQuadrature(std::shared_ptr<Polynomial> poly_in)
                                : poly(poly_in) {}

void GaussQuadrature::Calculate() {

  // Where do we get polyOrder?
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(polyOrder,polyOrder);

  for (unsigned int i=1; i<polyOrder+1; i++) {

    // Need to connect to the a,b, and c functions for the specic polynomial
    double alpha_i = -b(i)/c(i);
    double beta_i = std::sqrt(c(i+1)/(a(i)*a(i+1)));

    // Diagonal entry of J
    J(i,i) = alpha_i;

    // Off diagonal entries of J
    if (i < polyOrder) {
      J(i,i+1) = beta_i;
      J(i+1,i) = beta_i;
    }

  }

  EigenSolver<MatrixXd> es(J, true);

  // Set gauss points
  // Will the eigenvalues function return an Eigen::VectorXd?
  gauss_pts = es.eigenvalues();

  // Get mu_0 value (integral of weighting function)
  // Where do we get this from?
  double mu_0 = 1;

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
