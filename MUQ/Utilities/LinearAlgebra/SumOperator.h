#ifndef SUMOPERATOR_H
#define SUMOPERATOR_H

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"


namespace muq
{
namespace Utilities
{

    class SumOperator : public LinearOperator
    {
    public:
        SumOperator(std::shared_ptr<LinearOperator> Ain,
                    std::shared_ptr<LinearOperator> Bin);

        /** Apply the linear operator to a vector */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Apply the transpose of the linear operator to a vector. */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

    private:
        std::shared_ptr<LinearOperator> A, B;

    }; // class SumOperator


} // namespace Utilities
} // namespace muq

#endif // #ifndef SUMOPERATOR_H
