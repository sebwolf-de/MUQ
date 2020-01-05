#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"

#include "MUQ/Modeling/LinearAlgebra/LOBPCG.h"

#include <random>

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MUQ/Utilities/Exceptions.h"
#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Modeling;
using namespace muq::Utilities;

TEST(LOBPCG, Diagonal_BigBlock)
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
    const double solveTol = 1e-7;
    const double eigTol = 0.0;
    const int blockSize = numEigs;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);

    solver.compute(op);

    EXPECT_NEAR(A(0,0), solver.eigenvalues()(0), 1e-4);
    EXPECT_NEAR(A(1,1), solver.eigenvalues()(1), 1e-4);
    EXPECT_NEAR(A(2,2), solver.eigenvalues()(2), 1e-4);
    EXPECT_NEAR(A(3,3), solver.eigenvalues()(3), 1e-4);
}

TEST(LOBPCG, Diagonal_MedBlock)
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
    const double solveTol = 1e-7;
    const double eigTol = 0.0;
    const int blockSize = 2;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);

    solver.compute(op);

    EXPECT_NEAR(A(0,0), solver.eigenvalues()(0), 1e-4);
    EXPECT_NEAR(A(1,1), solver.eigenvalues()(1), 1e-4);
    EXPECT_NEAR(A(2,2), solver.eigenvalues()(2), 1e-4);
    EXPECT_NEAR(A(3,3), solver.eigenvalues()(3), 1e-4);
}

TEST(LOBPCG, Diagonal_SmallBlock)
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
    const double solveTol = 1e-7;
    const double eigTol = 0.0;
    const int blockSize = 1;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);

    solver.compute(op);

    EXPECT_NEAR(A(0,0), solver.eigenvalues()(0), 1e-4);
    EXPECT_NEAR(A(1,1), solver.eigenvalues()(1), 1e-4);
    EXPECT_NEAR(A(2,2), solver.eigenvalues()(2), 1e-4);
    EXPECT_NEAR(A(3,3), solver.eigenvalues()(3), 1e-4);
}

TEST(LOBPCG, Diagonal_Tolerance)
{

    const int dim = 20;
    const double nugget = 1e-12;

    // Create a random symmetric positive definite matrix
    Eigen::MatrixXd A = nugget*Eigen::MatrixXd::Identity(dim,dim);
    A(0,0) = 1.0;
    A(1,1) = 0.5;
    A(2,2) = 0.09;
    A(3,3) = 0.05;

    auto op = LinearOperator::Create(A);

    const int numEigs = 4;
    const double solveTol = 1e-7;
    double relTol = 0.1;
    double absTol = 0.0;
    const int blockSize =1;

    {
      LOBPCG solver(numEigs, relTol, absTol, blockSize, solveTol);

      solver.compute(op);

      EXPECT_EQ(3, solver.eigenvalues().size());
      EXPECT_NEAR(A(0,0), solver.eigenvalues()(0), 1e-4);
      EXPECT_NEAR(A(1,1), solver.eigenvalues()(1), 1e-4);
      EXPECT_NEAR(A(2,2), solver.eigenvalues()(2), 1e-4);
    }

    {
      double relTol = 0.0;
      double absTol = 0.1;

      LOBPCG solver(numEigs, relTol, absTol, blockSize, solveTol);

      solver.compute(op);

      EXPECT_EQ(3, solver.eigenvalues().size());
      EXPECT_NEAR(A(0,0), solver.eigenvalues()(0), 1e-4);
      EXPECT_NEAR(A(1,1), solver.eigenvalues()(1), 1e-4);
      EXPECT_NEAR(A(2,2), solver.eigenvalues()(2), 1e-4);
    }
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
    const double solveTol = 1e-7;
    const double eigTol = 0.0;
    const int blockSize = numEigs;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);
    solver.compute(op);

    for(unsigned int i=0; i<numEigs; ++i)
      EXPECT_NEAR(trueVals(subDim-1-i), solver.eigenvalues()(i), 1e-10);

}

TEST(LOBPCG, RandomIdentity)
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


    auto opA = LinearOperator::Create(A);
    auto opB = std::make_shared<IdentityOperator>(dim);

    const int numEigs = subDim-2;
    const double solveTol = 1e-7;
    const double eigTol = 0.0;
    const int blockSize = numEigs;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);

    solver.compute(opA,opB);

    for(unsigned int i=0; i<numEigs; ++i)
      EXPECT_NEAR(trueVals(subDim-1-i), solver.eigenvalues()(i), 1e-10);

}

TEST(LOBPCG, RandomGeneral)
{
    //RandomGenerator::SetSeed(2012);

    // The matrix size
    const int dim = 20;
    const int subDim = 3;
    const double nugget = 1e-12;

    // Generate two random positive definite matrices to define the Generalized eigenvalue problem
    Eigen::MatrixXd temp = Eigen::MatrixXd::Random(dim,subDim);
    Eigen::MatrixXd A = temp * temp.transpose() + nugget*Eigen::MatrixXd::Identity(dim,dim);

    temp = Eigen::MatrixXd::Random(dim,10);
    Eigen::MatrixXd B = temp * temp.transpose() + 1e-1*Eigen::MatrixXd::Identity(dim,dim);

    // Solve the generalized problem with Eigen3's direct solver for comparison
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> dirSolver(A,B);

    // Now compute the solution with MUQ's LOBPCG implementation
    auto opA = LinearOperator::Create(A);
    auto opB = LinearOperator::Create(B);

    const int numEigs = 3;
    const double solveTol = 1e-6;
    const double eigTol = 0.0;
    const int blockSize = numEigs;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);
    solver.compute(opA,opB);

    for(int i=1;i<numEigs+1;++i)
      EXPECT_NEAR(dirSolver.eigenvalues()(dim-i), solver.eigenvalues()(i-1),1e-3*dirSolver.eigenvalues()(dim-i));
}

TEST(LOBPCG, RandomGeneralPrecond)
{
    //RandomGenerator::SetSeed(2012);

    // The matrix size
    const int dim = 100;
    const int subDim = 5;
    const double nugget = 1e-12;

    // Generate two random positive definite matrices to define the Generalized eigenvalue problem
    Eigen::MatrixXd temp = Eigen::MatrixXd::Random(dim,subDim);
    Eigen::MatrixXd A = temp * temp.transpose() + nugget*Eigen::MatrixXd::Identity(dim,dim);

    temp = Eigen::MatrixXd::Random(dim,10);
    Eigen::MatrixXd B = temp * temp.transpose() + 1e-3*Eigen::MatrixXd::Identity(dim,dim);

    // Solve the generalized problem with Eigen3's direct solver for comparison
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> dirSolver(A,B);

    // Now compute the solution with MUQ's LOBPCG implementation
    auto opA = LinearOperator::Create(A);
    auto opB = LinearOperator::Create(B);
    auto opM = LinearOperator::Create(B.ldlt().solve(Eigen::MatrixXd::Identity(dim,dim)).eval());

    const int numEigs = 5;
    const double solveTol = 1e-6;
    const double eigTol = 0.0;
    const int blockSize = numEigs;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);

    solver.compute(opA,opB, opM);

    for(int i=0;i<numEigs;++i)
      EXPECT_NEAR(dirSolver.eigenvalues()(dim-i-1), solver.eigenvalues()(i),1e-3*dirSolver.eigenvalues()(dim-i-1));
}

TEST(LOBPCG, RandomGeneralPrecondConstrained)
{
    //RandomGenerator::SetSeed(2012);

    // The matrix size
    const int dim = 100;
    const int subDim = 5;
    const double nugget = 1e-12;

    // Generate two random positive definite matrices to define the Generalized eigenvalue problem
    Eigen::MatrixXd temp = Eigen::MatrixXd::Random(dim,subDim);
    Eigen::MatrixXd A = temp * temp.transpose() + nugget*Eigen::MatrixXd::Identity(dim,dim);

    temp = Eigen::MatrixXd::Random(dim,10);
    Eigen::MatrixXd B = temp * temp.transpose() + 1e-3*Eigen::MatrixXd::Identity(dim,dim);

    // Solve the generalized problem with Eigen3's direct solver for comparison
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> dirSolver(A,B);

    // Now compute the solution with MUQ's LOBPCG implementation
    auto opA = LinearOperator::Create(A);
    auto opB = LinearOperator::Create(B);
    auto opM = LinearOperator::Create(B.ldlt().solve(Eigen::MatrixXd::Identity(dim,dim)).eval());

    Eigen::MatrixXd constraints = dirSolver.eigenvectors().rightCols(2);


    const int numEigs = 3;
    const double solveTol = 1e-6;
    const double eigTol = 0.0;
    const int blockSize = numEigs;

    LOBPCG solver(numEigs, eigTol, eigTol, blockSize, solveTol);

    solver.compute(opA, constraints, opB, opM);

    for(int i=0;i<numEigs;++i)
      EXPECT_NEAR(dirSolver.eigenvalues()(dim-i-constraints.cols()-1), solver.eigenvalues()(i),1e-3*dirSolver.eigenvalues()(dim-i-constraints.cols()-1));
}
