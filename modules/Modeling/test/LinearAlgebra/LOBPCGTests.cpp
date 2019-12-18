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


TEST(LOBPCG, Random)
{
    // The matrix size
    const int dim = 30;

    // The dimension of the range of the matrix
    const int subDim = 5;

    // The true nonzero eigenvalues
    Eigen::VectorXd trueVals = Eigen::VectorXd::Random(subDim).array().abs();
    std::sort(trueVals.data(), trueVals.data()+subDim);

    // The true eigenvectors
    Eigen::MatrixXd trueVecs = Eigen::MatrixXd::Random(dim,subDim);

    // Orthonormalize the vectors
    for(unsigned int i=0; i<subDim; ++i){
      trueVecs.col(i) /= trueVecs.col(i).norm();

      for(unsigned int j=i+1; j<subDim; ++j){
        trueVecs.col(j) -= trueVecs.col(i).dot(trueVecs.col(j))*trueVecs.col(i);
      }
    }

    const double nugget = 1e-12;

    // Create a matrix with the specified eigenvalues and eigenvectors
    Eigen::MatrixXd A = trueVecs * trueVals.asDiagonal() * trueVecs.transpose() + nugget*Eigen::MatrixXd::Identity(dim,dim);


    auto op = LinearOperator::Create(A);

    const int numEigs = subDim-2;
    const double tol = 1e-7;
    LOBPCG solver(numEigs, tol);

    solver.compute(op);

    for(unsigned int i=0; i<subDim-2; ++i)
      EXPECT_NEAR(trueVals(subDim-i-1), solver.eigenvalues()(subDim-3-i), 1e-10);

}
