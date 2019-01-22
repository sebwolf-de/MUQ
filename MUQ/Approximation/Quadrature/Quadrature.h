#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <memory>

#include <Eigen/Core>

namespace muq {
  namespace Approximation {

    /**
    @defgroup Quadrature
    @brief Tools for constructing multivariate quadrature rules.
    @details A whole bunch of details and examples to come...
    */

    /**
     @class Quadrature
     @ingroup Quadrature
     @brief Base class for multivariate quadrature rules.
     @detail An abstract class for computing nodes and weights of general quadrature rules.
     @seealso GaussQuadrature
     */
    class Quadrature {
    public:

      Quadrature(unsigned int dimIn) : dim(dimIn){};


      virtual ~Quadrature() = default;

      virtual void Compute(unsigned int order) = 0;

      /** Base implementation of Compute.  Assumes the quadrature rule is 1d
          and then calls Compute with the first component of the orders vector.

          Multivariate quadrature rules should override this function.
      */
      virtual void Compute(Eigen::RowVectorXi const& orders) { assert(orders.size()==1); Compute(orders(0));}

      /** Return the dimension of the quadrature rule. */
      virtual unsigned int Dim() const{return dim;};

      /** Return the quadrature points.  The output is (Dim x NumPts), where
          dim is the dimension of the integral under consideration and NumPts
          is the number of quadrature points used in the rule.
      */
      virtual Eigen::MatrixXd const& Points() const{return pts;};

      /** Return the quadrature weights.
      */
      virtual Eigen::VectorXd const& Weights() const{return wts;};


    protected:
      unsigned int dim;
      Eigen::MatrixXd pts;
      Eigen::VectorXd wts;
    };
  }
}

#endif /* QUADRATURE_H_ */
