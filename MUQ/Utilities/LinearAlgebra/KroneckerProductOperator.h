#ifndef KRONECKERPRODUCTOPERATOR_H
#define KRONECKERPRODUCTOPERATOR_H

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"

namespace muq
{
namespace Utilities
{

    /** @class KroneckerProductOperator
        @ingroup LinearOperators
        @brief Defines the Kronecker product of two other linear operators.
        @details Let \f$A\in\mathbb{R}^M\times\mathbb{R}^N\f$ and \f$B\f$ be two rectangular matrices.  The Kronecker product of \f$A\f$ and \f$B\f$ is given in block form by
\f[
A\times B = \left[ \begin{array}{ccc} A_{11} B & \cdots & A_{1N} B\\ \vdots & \ddots & \vdots \\ A_{M1} & \cdots & A_{MN} B \end{array}\right],
\f]
    */
    class KroneckerProductOperator : public LinearOperator
    {

    public:
        KroneckerProductOperator(std::shared_ptr<LinearOperator> Ain,
                                 std::shared_ptr<LinearOperator> Bin);


        /** Fills in the reference \f$y\f$ with \f$y=Ax\f$ */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Fill in the reference \f$y\f$ with \f$y = A^Txf$ */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

    private:
        std::shared_ptr<LinearOperator> A;
        std::shared_ptr<LinearOperator> B;

    }; // class KroneckerProductOperator


} // namespace Utilities
} // namespace muq


#endif // #ifndef KRONECKERPRODUCTOPERATOR_H
