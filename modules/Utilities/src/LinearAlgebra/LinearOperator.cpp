#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"


using namespace muq::Utilities;

void LinearOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y)
{
    assert(y.cols()==x.cols());
    y = Apply(x);
};


void LinearOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y)
{
    assert(y.cols()==x.cols());
    y = ApplyTranspose(x);
};

Eigen::MatrixXd LinearOperator::GetMatrix()
{

    Eigen::MatrixXd output(nrows, ncols);
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Identity(ncols, ncols);

    for(int i=0; i<ncols; ++i)
        output.col(i) = Apply(rhs.col(i));

    return output;
}
