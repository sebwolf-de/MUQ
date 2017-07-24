#include "MUQ/Utilities/LinearAlgebra/CompanionMatrix.h"

using namespace muq::Utilities;



Eigen::MatrixXd CompanionMatrix::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    assert(x.rows() == ncols);
    
    Eigen::MatrixXd output(nrows, x.cols());
    
    output.block(0, 0, nrows-1, x.cols()) = x.block(1,0,nrows-1, x.cols());
    output.row(nrows-1) = lastRow.transpose()*x;

    return output;
}


Eigen::MatrixXd CompanionMatrix::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    assert(x.rows() == nrows);
    
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(ncols, x.cols());

    output.block(1, 0, ncols-1, x.cols()) = x.block(0,0,ncols-1,x.cols());
    output += lastRow * x.row(x.rows()-1);
    
    return output;
}
