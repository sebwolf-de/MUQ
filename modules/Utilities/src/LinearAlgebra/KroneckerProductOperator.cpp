#include "MUQ/Utilities/LinearAlgebra/KroneckerProductOperator.h"


using namespace muq::Utilities;

KroneckerProductOperator::KroneckerProductOperator(std::shared_ptr<LinearOperator> Ain,
                                                   std::shared_ptr<LinearOperator> Bin) : LinearOperator(Ain->rows()*Bin->rows(), Ain->cols()*Bin->cols()), A(Ain), B(Bin)
{

}


Eigen::MatrixXd KroneckerProductOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{

    Eigen::MatrixXd output(nrows, x.cols());
    
    for(int i=0; i<x.cols(); ++i)
    {
        Eigen::VectorXd xVec = x.col(i);
        Eigen::Map<Eigen::MatrixXd> xMat(xVec.data(), B->cols(), A->cols());
        Eigen::Map<Eigen::MatrixXd> bMat(&output(0,i), B->rows(), A->rows());

        bMat = A->Apply( B->Apply( xMat ).transpose() ).transpose();
    }

    return output;
}



Eigen::MatrixXd KroneckerProductOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    Eigen::MatrixXd output(ncols, x.cols());
    
    for(int i=0; i<x.cols(); ++i)
    {
        Eigen::VectorXd xVec = x.col(i);
        Eigen::Map<const Eigen::MatrixXd> xMat(xVec.data(), B->rows(), A->rows());
        Eigen::Map<Eigen::MatrixXd> bMat(&output(0,i), B->cols(), A->cols());

        bMat = A->ApplyTranspose( B->ApplyTranspose( xMat ).transpose() ).transpose();
    }

    return output;
}
