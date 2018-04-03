#ifndef GAUSSQUADRATURE_H_
#define GAUSSQUADRATURE_H_

// What header files do we need? Add here.

namespace muq {

  namespace Approximation {

    class GaussQuadrature {

    public:

      GaussQuadrature();

      GaussQuadrature(std::shared_ptr<Polynomial> poly_in);

      // Do we need a destructor?

      void Calculate();

    private:

      int polyOrder;

      Eigen::VectorXd gauss_pts;

      Eigen::VectorXd gauss_wts;

      std::shared_ptr<Polynomial> poly;

    };

  }

}

#endif
