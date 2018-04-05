#ifndef GAUSSQUADRATURE_H_
#define GAUSSQUADRATURE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace muq {

  namespace Approximation {

    class GaussQuadrature {

    public:

      GaussQuadrature();

      GaussQuadrature(std::shared_ptr<OrthogonalPolynomial> polyIn,
                      int polyOrderIn);

      void Calculate();

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
