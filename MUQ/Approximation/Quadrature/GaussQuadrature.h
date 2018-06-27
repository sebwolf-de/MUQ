#ifndef GAUSSQUADRATURE_H_
#define GAUSSQUADRATURE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace muq {

  namespace Approximation {

    /** @ingroup Polynomials
        @brief Class for computing Gauss Quadrature rules from an orthogonal polynomial family.
        @details Uses the Golub-Welsch algorithm to construct a Gauss Quadrature rule of a
        specified order.
    */
    class GaussQuadrature {

    public:

      GaussQuadrature();

      GaussQuadrature(std::shared_ptr<OrthogonalPolynomial> polyIn,
                      int polyOrderIn);

      void Compute();

      Eigen::VectorXd const& Points() const;

      Eigen::VectorXd const& Weights() const;

    private:

      Eigen::VectorXd gaussPts;

      Eigen::VectorXd gaussWts;

      std::shared_ptr<OrthogonalPolynomial> poly;

      int polyOrder;

    };

  }

}

#endif
