#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <memory>

#include <Eigen/Core>

namespace muq {
  namespace Approximation {

    /**
    @defgroup Quadrature
    @brief Tools for constructing univariate and multivariate quadrature rules.
    @details
    <b> Gauss Quadrature </b>
    One dimensional quadrature rules attempt to approximate a weighted integral
    of the form
    \f[
    \int_\Omega f(x) w(x) dx,
    \f]
    where \f$f(x):\mathbb{R}\rightarrow\mathbb{R}\f$ is the function we wish to
    integrate, \f$w(x):\mathbb{R}\rightarrow\mathbb{R}\f$ is a given
    weight function, and \f$\Omega\f$ defines the support of the weight function.
    To approximate this type of integral, we seek an approximation of the form
    \f[
    \int_\Omega f(x) w(x) dx \approx \sum_{i=1}^N w_i f(x_i),
    \f]
    where \f$w_i\f$ and \f$x_i\f$ are called the quadrature weights and quadrature
    points.  Together these are called the quadrature rule.

    Gauss quadrature refers to a particular way of defining quadrature rules
    that can exactly integrate polynomials of order \f$2N-1\f$.  The \f$N\f$ point
    Gauss quadrature rule for a weight function \f$w(x)\f$ corresponds to the roots
    of the degree-\f$N\f$ polynomial in a polynomial family that is orthogonal
    with respect o the weight function \f$w(x)\f$.  For example, Legendre polynomials
    are orthogonal with respect to \f$w(x) = I([-1,1])\f$, where \f$I([-1,1])\f$
    is the indicator function for the set \f$[-1,1]\f$.  The \f$N\f$ point
    Gauss-Legendre rule will then use the roots of the \f$N^{th}\f$ order Legendre
    polynomial to exactly integrate
    \f[
    \int_{-1}^1 f(x) dx,
    \f]
    when \f$f(x)\f$ is a polynomial with degree \f$\leq N\f$.  The following table
    matches orthogonal polynomial families with their corresponding weight function
    and MUQ class name.

    Weight Function \f$w(x)\f$  | Domain       | Polynomial Family | MUQ Class
    --------------------------- | ------------ | ------------------| ---------------
    \f$1\f$                     | \f$[-1,1]\f$ | Legendre          | `muq::Approximation::Legendre`
    \f$\exp\left[-\frac{1}{2}x^2\right] \f$ | \f$(-\infty, \infty)\f$ | Probabilist Hermite | `muq::Approximation::ProbabilistHermite`
    \f$\exp\left[-x^2\right] \f$ | \f$(-\infty, \infty)\f$ | Physicist Hermite | `muq::Approximation::PhysicistHermite`
    \f$\exp\left[-x\right] \f$  | \f$[0, \infty)\f$ | Laguerre | `muq::Approximation::Laguerre`
    \f$(1-x)^a (1+x)^b          | \f$[-1,1]\f$      | Jacobi   | `muq::Approximation::Jacobi`

    To construct a Gauss quadrature rule, MUQ requires the polynomial family to be
    defined first.  See below for examples.  Note that behind the scenes, MUQ
    uses the Golub-Welsch algorithm (see "Calculation of Quadrature Rules" Golub and Welsch 1969),
    to compute Gauss quadrature rules.

    Example of constructing a Gauss-Legendre rule:
    @code{.cpp}
    auto poly = std::make_shared<Legendre>();
    GaussQuadrature gaussLegendre(poly);

    gaussLegendre.Compute(order);
    Eigen::VectorXd gaussPts = gq.Points().transpose();
    Eigen::VectorXd gaussWts = gq.Weights();
    @endcode

    Set up for Gauss-Hermite rules:
    @code{.cpp}
    auto poly = std::make_shared<PhysicistHermite>();
    GaussQuadrature gaussHermite(poly);
    @endcode
    @code{.cpp}
    auto poly = std::make_shared<ProbabilistHermite>();
    GaussQuadrature gaussHermite(poly);
    @endcode

    Set up for Gauss-Laguerre
    @code{.cpp}
    auto poly = std::make_shared<Laguerre>();
    GaussQuadrature gaussLaguerre(poly);
    @endcode

    Set up for Gauss-Jacobi
    @code{.cpp}
    double a = 0.5;
    double b = 0.5;
    auto poly = std::make_shared<Jacobi>(a,b);
    GaussQuadrature gaussJacobi(poly);
    @endcode

    Note that it is also possible to choose a polynomial family by passing
    the polynomial name as a string to the `OrthogonalPolynomial::Construct`
    method.  For example,
    @code{.cpp}
    std::string polyName = "Legendre";
    auto poly = OrthogonalPolynomial::Construct(polyName);
    GaussQuadrature gaussQuad(poly);
    @endcode

    <b> Other 1D Quadrature Rules </b>
    The nexted Clenshaw Curtis quadrature:
    @code{.cpp}
    unsigned int numPts = 10;
    ClenshawCurtisQuadrature ccQuad;
    ccQuad.Compute(numPts);
    @endcode
    
    <b> Tensor Product Rules </b>

    <b> Smolyak Rules </b>
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
