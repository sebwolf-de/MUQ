#ifndef MONOTONEEXPANSION_H
#define MONOTONEEXPANSION_H

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include <boost/property_tree/ptree.hpp>

namespace muq{
  namespace Approximation{

    /** @class MonotoneExpansion
        @ingroup Polynomials
        @brief Defines a parameterized family of multivariate lower triangular monotone functions.
        @details This class defines functions that take the form
        \f[ f_i(x) = g_i(x_1, \ldots, x_{i-1}) + \int_{-\infty}^{x_i} \left[ h_i(x_1,\ldots, x_{i-1}, y)\right]^2 + \epsilon dy, \f]
        where the functions \f$g_i\f$ and \f$h_i\f$ are both represented through the
        muq::Approximation::BasisExpansion class, and \f$\epsilon>0\f$ is a small nugget
        ensuring that \f$f_i\f$ is strictly increasing with \f$x_i\f$.  By squaring \f$h_i\f$,
        we guarantee that the integrand is positive, which ensures that the integral itself increases
        with \f$x_i\f$, i.e., \f$\partial f_i / \partial x_i > 0 \f$.  Also notice how
        the \f$i^{th}\f$ output of the function \f$f\f$ only depends on the first \f$i\$ inputs,
        so the function \f$f\f$ is lower triangular. Combined with the fact that \f$\partial f_i / \partial x_i > 0 \f$,
        this means that \f$f\f$ will be invertible.
    */
    class MonotoneExpansion : public muq::Modeling::WorkPiece{

    public:
      //MonotoneExpansion(boost::property_tree::ptree & params);

      MonotoneExpansion(std::shared_ptr<BasisExpansion> monotonePartsIn);

      MonotoneExpansion(std::vector<std::shared_ptr<BasisExpansion>> const& generalPartsIn,
                        std::vector<std::shared_ptr<BasisExpansion>> const& monotonePartsIn);

      /** Returns the number of coefficients across all expansions in all dimension, i.e.,
          the total number of degrees of freedom describing this expansion.
      */
      unsigned NumTerms() const;

      /** Get the current expansion coefficients.  The coefficients are ordered
          with coefficients from the general parts preceding coefficients for
          the monotone parts.  Within the general and montone parts, the coefficients
          are ordered according to the output dimension.
      */
      virtual Eigen::VectorXd GetCoeffs() const;

      virtual void SetCoeffs(Eigen::VectorXd const& allCoeffs);

      /** Returns the log determinant of the Jacobian matrix (wrt x) at a particular point. */
      virtual double LogDeterminant(Eigen::VectorXd const& evalPt);
      virtual double LogDeterminant(Eigen::VectorXd const& evalPt, Eigen::VectorXd const& coeffs);

      /** Returns the gradient of the Jacobian log determinant with respect to the coefficients. */
      virtual Eigen::VectorXd GradLogDeterminant(Eigen::VectorXd const& evalPt);
      virtual Eigen::VectorXd GradLogDeterminant(Eigen::VectorXd const& evalPt,
                                                 Eigen::VectorXd const& coeffs);

    protected:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      virtual void JacobianImpl(unsigned int const                           wrtIn,
                                unsigned int const                           wrtOut,
                                muq::Modeling::ref_vector<boost::any> const& inputs) override;

      std::vector<std::shared_ptr<BasisExpansion>> generalParts;
      std::vector<std::shared_ptr<BasisExpansion>> monotoneParts;

      Eigen::VectorXd quadWeights;
      Eigen::VectorXd quadPts; // quadrature points on [0,1]

    }; // class MonotoneExpansion

  } // namespace Approximation
} // namespace muq


#endif
