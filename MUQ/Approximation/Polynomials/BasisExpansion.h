#ifndef BASISEXPANSION_H
#define BASISEXPANSION_H

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

namespace muq{
  namespace Approximation{

    class MonotoneExpansion;
    
    /** @class BasisExpansion
        @ingroup Polynomials
        @brief Class for defining expansions of basis functions defined by a
        MultiIndexSet and collection of IndexScalarBasis functions.
        @details Consider an expansion of the form
                 \f[f_j(x) = \sum_{\alpha\in A} c_{j,\alpha} \Phi_\alpha(x), \f]
                 where \f$x\in \mathbb{R}^N\f$, \f$c_{j,\alpha} \in \mathbb{R}\f$,
                 is a coefficient, \f$\alpha\f$ is a multindex in the set of
                 indices \f$A\f$, and \f$\Phi_\alpha(x)\f$ is a multivariate
                 basis function defined by the multiindex, which takes the form
                 \f[ \Phi_\alpha(x) = \Prod_{i=1}^N \phi_{i}(x_i,\alpha_i). \f]
                 The univariate functions \f$\phi_{i}(x_i,\alpha_i)\f$ can be
                 polynomials, Hermite functions, or some ther IndexScalarBasis.
                 For example, we could use Hermite polynomials for \f$i=0\f$ and
                 Legendre polynomials for \f$i=1\f$, so that \f$\phi_0(x_0,\alpha_0)\f$
                 would be a Hermite polynomial of order \f$\alpha_0\f$ and
                 \f$\phi_1(x_1,\alpha_1)\f$ would be Legendre polynomial of order
                 \f$\alpha_1\f$.

                 Evaluating this WorkPiece should be done through the muq::Modeling::WorkGraph
                 interface, (WorkPiece::Evaluate, WorkPiece::Jacobian, etc...).  The input
                 arguments is either the vector \f$x\f$ or both the vector \f$x\f$ and a MatrixXd
                 of coefficients defining \f$c_{j,\alpha}\f$ in the expansion.   If the coefficients
                 are not passed, the most recently set coefficients are used.  For example, `output1`
                 and `output2` will be the same below:
                 @code
auto expansion = std::make_shared<BasisExpansion>();

Eigen::VectorXd x;
Eigen::MatrixXd c;

// Fill in x and c ...

boost::any output1 = expansion->Evaluate(x,c)[0];
Eigen::MatrixXd outputVec1 = boost::any_cast<Eigen::MatrixXd>(output1);

boost::any output2 = expansion->Evaluate(x)[0];

Eigen::MatrixXd outputVec2 = boost::any_cast<Eigen::MatrixXd>(output1);

                 @endcode
        @seealso muq::Utilities::MultiIndexSet, muq::Utilities::IndexScalarBasis, muq::Modeling::WorkPiece
    */
    class BasisExpansion : public muq::Modeling::WorkPiece{

      friend class MonotoneExpansion;

    public:

      /** Construct expansion by specifying basis components (i.e., the family for each \f$\phi_i\f$).
          Initializes the expansion to a single constant (i.e., \f$\alpha=0\f$) term with coefficient \f$0\f$.
      */
      BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn);

      /** Construct the expansion by specifying both the basis families and multi-indices.
          Sets all coefficients in the expansion to zero.
      */
      BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                     std::shared_ptr<MultiIndexSet>                          multisIn);

      /** Construct the expansion by specifying all ingredients: the basis family, multi-indices, and coefficients. */
      BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                     std::shared_ptr<MultiIndexSet>                          multisIn,
                     Eigen::MatrixXd                                  const& coeffsIn);

      virtual ~BasisExpansion() = default;


      /** Get the \f$k^{th}\f$ derivative of WorkPiece with respect to the input \f$x\f$. */
      virtual Eigen::MatrixXd Derivative(unsigned                            derivOrder,
                                         std::vector<Eigen::MatrixXd> const& inputs);

    protected:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      virtual void JacobianImpl(unsigned int const                           wrtIn,
                                unsigned int const                           wrtOut,
                                muq::Modeling::ref_vector<boost::any> const& inputs) override;

      // Evaluates all the terms in the expansion, but does not multiply by coefficients
      Eigen::VectorXd GetAllTerms();

      // MultiIndexSet defining each term in the expansion
      std::shared_ptr<MultiIndexSet> multis;

      // Coefficients for the output
      Eigen::MatrixXd coeffs;

      // Components of the basis functions
      std::vector<std::shared_ptr<IndexedScalarBasis>> basisComps;

    };

  }
}

#endif
