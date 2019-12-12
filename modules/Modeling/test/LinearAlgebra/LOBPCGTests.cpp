#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"

#include "MUQ/Modeling/LinearAlgebra/LOBPCG.h"

#include <random>

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

TEST(LOBPCG, Diagonal)
{

    const int dim = 20;
    const double nugget = 1e-12;

    // Create a random symmetric positive definite matrix
    Eigen::MatrixXd A = nugget*Eigen::MatrixXd::Identity(dim,dim);
    A(0,0) = 1.0;
    A(1,1) = 0.5;
    A(2,2) = 0.25;
    A(3,3) = 0.125;

    auto op = LinearOperator::Create(A);

    const int numEigs = 4;
    const double tol = 1e-7;
    LOBPCG solver(numEigs, tol);

    solver.compute(op);
}
