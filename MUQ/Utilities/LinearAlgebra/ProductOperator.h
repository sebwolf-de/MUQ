#ifndef PRODUCTOPERATOR_H
#define PRODUCTOPERATOR_H

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"


namespace muq
{
namespace Utilities
{

    /** @class ProductOperator
        @ingroup LinearOperators
        @brief The product of two linear operators, \f$C=A*B\f$
    */
    class ProductOperator : public LinearOperator
    {
    public:
        ProductOperator(std::shared_ptr<LinearOperator> Ain,
                        std::shared_ptr<LinearOperator> Bin);

        virtual ~ProductOperator(){};
        
        /** Apply the linear operator to a vector */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Apply the transpose of the linear operator to a vector. */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

        virtual Eigen::MatrixXd GetMatrix() override;
        
    private:
        std::shared_ptr<LinearOperator> A, B;

    }; // class SumOperator


} // namespace Utilities
} // namespace muq

#endif // #ifndef DIAGONALOPERATOR_H
