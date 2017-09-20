#ifndef DIAGONALOPERATOR_H
#define DIAGONALOPERATOR_H

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"


namespace muq
{
namespace Utilities
{

    class DiagonalOperator : public LinearOperator
    {
    public:
        DiagonalOperator(Eigen::VectorXd const& diagIn);

        virtual ~DiagonalOperator(){};
        
        /** Apply the linear operator to a vector */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Apply the transpose of the linear operator to a vector. */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override{return Apply(x);};

        virtual Eigen::MatrixXd GetMatrix() override{return diag.asDiagonal();};
        
    private:
        const Eigen::VectorXd diag;

    }; // class SumOperator


} // namespace Utilities
} // namespace muq

#endif // #ifndef DIAGONALOPERATOR_H
