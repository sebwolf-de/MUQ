#include "MUQ/Utilities/LinearAlgebra/DiagonalOperator.h"

using namespace muq::Utilities;


DiagonalOperator::DiagonalOperator(Eigen::VectorXd const& diagIn) : LinearOperator(diagIn.rows(), diagIn.rows())
{}

/** Apply the linear operator to a vector */
Eigen::MatrixXd DiagonalOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    return diag.asDiagonal()*x;
}
