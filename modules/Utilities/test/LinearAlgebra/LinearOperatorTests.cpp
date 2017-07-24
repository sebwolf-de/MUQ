#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"
#include "MUQ/Utilities/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Utilities/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"

#include <random>

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace muq::Utilities;

TEST(Utilties_LinearOperator, EigenDense)
{

    const int rows = 5;
    const int cols = 2;
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(rows,cols);

    auto denseOp = LinearOperator::Create(A);

    Eigen::VectorXd x = Eigen::VectorXd::Ones(cols);
    Eigen::VectorXd ytrue, ytest;

    ytrue = A*x;
    ytest = denseOp->Apply(x);

    for(int i=0; i<rows; ++i)
        EXPECT_DOUBLE_EQ(ytrue(i), ytest(i));

    // TRANSPOSE TESTS
    x = Eigen::VectorXd::Ones(rows);
    ytrue = A.transpose()*x;
    ytest = denseOp->ApplyTranspose(x);

    for(int i=0; i<cols; ++i)
        EXPECT_DOUBLE_EQ(ytrue(i), ytest(i));

}

TEST(Utilties_LinearOperator, EigenSparse)
{

    const int nnz = 10;
    const int rows = 20;
    const int cols = 10;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> rowRG(0, rows-1);
    std::uniform_int_distribution<> colRG(0, cols-1);
    std::normal_distribution<> valRG(0.0,1.0);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(nnz);

    // Create a random sparse matrix with at most nnz non zero entries
    for(int i=0; i<nnz; ++i){
        
        int    randRow = rowRG(gen);
        int    randCol = colRG(gen);
        double randVal = valRG(gen);
        
        tripletList.push_back(T(randRow, randCol, randVal));
    }

    Eigen::SparseMatrix<double> A(rows,cols);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    
    auto denseOp = LinearOperator::Create(A);

    Eigen::VectorXd x = Eigen::VectorXd::Ones(cols);
    Eigen::VectorXd ytrue, ytest;

    ytrue = A*x;
    ytest = denseOp->Apply(x);

    for(int i=0; i<rows; ++i)
        EXPECT_DOUBLE_EQ(ytrue(i), ytest(i));

    // TRANSPOSE TESTS
    x = Eigen::VectorXd::Ones(rows);
    ytrue = A.transpose()*x;
    ytest = denseOp->ApplyTranspose(x);
    
    for(int i=0; i<cols; ++i)
        EXPECT_DOUBLE_EQ(ytrue(i), ytest(i));

}


TEST(Utilities_LinearOperator, Conversion)
{
    Eigen::MatrixXd A(10,4);
    std::shared_ptr<LinearOperator> Aop = LinearOperator::Create(A);

}


TEST(Utilities_LinearOperator, CompanionMatrix)
{

    Eigen::VectorXd lastRow(3);
    lastRow << 1.0, 1.5, 1.75;

    CompanionMatrix A(lastRow);
    
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(3,10);

    // Test the multiplication
    Eigen::MatrixXd b = A.Apply(x);

    Eigen::VectorXd trueRow = lastRow.transpose()*x;
        
    for(int j=0; j<x.cols(); ++j)
    {
        for(int i=0; i<b.rows()-1; ++i)
        {
            EXPECT_DOUBLE_EQ(x(i+1,j), b(i,j));
        }

        EXPECT_DOUBLE_EQ(trueRow(j), b(b.rows()-1,j));
    }

    // Test the transpose application
    b = A.ApplyTranspose(x);
    Eigen::MatrixXd trueTrans = Eigen::MatrixXd::Zero(3,3);
    trueTrans.block(1,0,2,2) = Eigen::MatrixXd::Identity(2,2);
    trueTrans.col(2) = lastRow;
    Eigen::MatrixXd trueB = trueTrans*x;
    
    for(int j=0; j<b.cols(); ++j)
    {
        for(int i=0; i<b.rows(); ++i)
        {
            EXPECT_DOUBLE_EQ(trueB(i,j), b(i,j));
        }
    }

}


TEST(Utilities_LinearOperator, BlockDiagonal)
{

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2,3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(4,2);

    std::vector<std::shared_ptr<LinearOperator>> ops(2);
    ops.at(0) = LinearOperator::Create(A);
    ops.at(1) = LinearOperator::Create(B);

    auto ABop = std::make_shared<BlockDiagonalOperator>(ops);

    Eigen::MatrixXd x = Eigen::MatrixXd::Random(A.cols()+B.cols(), 10);

    Eigen::MatrixXd bTrue(A.rows() + B.rows(), x.cols());
    bTrue.block(0,0, A.rows(), bTrue.cols()) = A * x.block(0, 0, A.cols(), x.cols());
    bTrue.block(A.rows(),0, B.rows(), bTrue.cols()) = B * x.block(A.cols(), 0, B.cols(), x.cols());
    

    Eigen::MatrixXd bOp = ABop->Apply(x);
    for(int j=0; j<bOp.cols(); ++j){
        for(int i=0; i<bOp.rows(); ++i)
            EXPECT_DOUBLE_EQ(bTrue(i,j), bOp(i,j));
    }


    x = Eigen::MatrixXd::Random(A.rows() + B.rows(), 10);
    bTrue.resize(A.cols() + B.cols(), x.cols());
    bTrue.block(0,0,A.cols(), bTrue.cols()) = A.transpose() * x.block(0,0,A.rows(), x.cols());
    bTrue.block(A.cols(),0,B.cols(), bTrue.cols()) = B.transpose() * x.block(A.rows(), 0 , B.rows(), x.cols());


    bOp = ABop->ApplyTranspose(x);
    for(int j=0; j<bOp.cols(); ++j){
        for(int i=0; i<bOp.rows(); ++i)
            EXPECT_DOUBLE_EQ(bTrue(i,j), bOp(i,j));
    }

}
