#ifndef IDENTITYOPERATOR_H_
#define IDENTITYOPERATOR_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    /**
      @brief A muq::Modeling::ModPiece that is equivalent to applying the identity matrix to a single input vector.
      @details This ModPiece simply return its input, but is useful for using a
      variable multiple times in a WorkGraph.  For example, for evaluate a function
      \f$f( g(x), h(x))\f$, you could define \f$x\f$ as an IdentityOperator node
      in a WorkGraph and then add edges from the output of \f$x\f$ to \f$g\f$ and \f$h\f$.

      In previous versions of MUQ, this was called the VectorPassthroughModel.
    */
    class IdentityOperator: public ModPiece {
    public:

      /**
	       @param[in] dim The dimension of the single input and output vector
       */
      IdentityOperator(unsigned int dim);


    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

    }; // class ConstantVector

  } // namespace Modeling
} // namespace muq

#endif
