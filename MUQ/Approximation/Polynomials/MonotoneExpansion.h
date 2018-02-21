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
        we guarantee that the integrand positive, which ensures that the integral itself increases
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


    protected:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      virtual void JacobianImpl(unsigned int const                           wrtIn,
                                unsigned int const                           wrtOut,
                                muq::Modeling::ref_vector<boost::any> const& inputs) override;

      std::vector<std::shared_ptr<BasisExpansion>> generalParts;
      std::vector<std::shared_ptr<BasisExpansion>> monotoneParts;

      Eigen::VectorXd coeffs;
      
      Eigen::VectorXd quadWeights;
      Eigen::VectorXd quadPts; // quadrature points on [0,1]

    }; // class MonotoneExpansion

  } // namespace Approximation
} // namespace muq


#endif
