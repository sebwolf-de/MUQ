#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;

TEST(Approximation_GP, VectorNorm)
{

    int dim = 100;
    Eigen::VectorXd v1 = Eigen::VectorXd::Random(dim);
    Eigen::VectorXd v2 = Eigen::VectorXd::Random(dim);

    EXPECT_DOUBLE_EQ((v2-v1).norm(), CalcDistance(v1,v2));
}

TEST(Approximation_GP, HyperFit1d)
{
    const unsigned numPred  = 50;
    const unsigned maxTrain = 50;
    
    const double pi = 4.0 * atan(1.0);

    // Set linearly space locations where we want to evaluate the Gaussian Process
    Eigen::RowVectorXd predLocs(numPred);
    predLocs.setLinSpaced(numPred, 0,1);

    // Generate random training locations
    Eigen::RowVectorXd trainLocs(maxTrain);
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    for(int i=0; i<maxTrain; ++i)
	trainLocs(i) = uniform_dist(e1);

    // Generate the training data (with some random noise)
    Eigen::RowVectorXd trainData(maxTrain);
    for(int i=0; i<maxTrain; ++i)
	trainData(i) = sin(4*2*pi*trainLocs(i) ) + sqrt(1e-4)*normal_dist(e1);

    const unsigned dim = 1;
    auto kernel = SquaredExpKernel(dim, 2.0, 0.35, {0.1,10} )*PeriodicKernel(dim, 1.0, 0.75, 0.25, {0.5,5.0}, {0.5,5.0}, {0.25,0.5}) + WhiteNoiseKernel(dim, 1e-3, {1e-8,100});    
   
    // Create the GP
    ConstantMean mean(dim, 1);
    auto gp = ConstructGP(mean, kernel);


    gp.Fit(trainLocs, trainData);
    
    // Make a prediction
    auto post = gp.Predict(predLocs);

    //for(int i=0; i<predLocs.cols(); ++i)
    //	std::cout << "[" << predLocs(0,i) << ", " << post.mean(0,i) << "]," << std::endl;
    
}

TEST(Approximation_GP, HyperFit2d)
{
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::normal_distribution<double> normal_dist(0.0, 1.0);


    const unsigned numPred  = 50;
    const unsigned maxTrain = 50;

    
    // Set linearly space locations where we want to evaluate the Gaussian Process
    Eigen::MatrixXd predLocs(2, numPred);
    for(int i=0; i<numPred; ++i)
    {
	predLocs(0,i) = uniform_dist(e1);
	predLocs(1,i) = uniform_dist(e1);
    }
    
    // Generate random training locations
    Eigen::MatrixXd trainLocs(2, maxTrain);
  
    for(int i=0; i<maxTrain; ++i)
    {
	trainLocs(0,i) = uniform_dist(e1);
	trainLocs(1,i) = uniform_dist(e1);
    }
    
    // Generate the training data (with some random noise)
    Eigen::RowVectorXd trainData(maxTrain);
    for(int i=0; i<maxTrain; ++i)
	trainData(i) = trainLocs.col(i).squaredNorm();

    // define a tensor product kernel
    std::vector<unsigned> inds1{0};
    std::vector<unsigned> inds2{1};
    const unsigned dim = 2;
    auto kernel = SquaredExpKernel(dim, inds1, 2.0, 0.35, {0.1,10} )*SquaredExpKernel(dim, inds2, 2.0, 0.35, {0.1,10} );

    
    // Create the GP
    ConstantMean mean(dim, 1);
    auto gp = ConstructGP(mean, kernel);

    gp.Fit(trainLocs, trainData);
    gp.Optimize();
    
    // Make a prediction
    auto post = gp.Predict(predLocs);

}

TEST(Approximation_GP, LinearTransformKernel)
{
    
    const unsigned dim = 2;
    Eigen::MatrixXd sigma2(2,2);
    sigma2 << 1.0, 0.9,
	      0.9, 1.5;
    
    auto kernel = ConstantKernel(dim, sigma2) * SquaredExpKernel(dim, 2.0, 0.35 );

    EXPECT_EQ(2, kernel.coDim);
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3,2);

    auto kernel2 = TransformKernel(A, kernel);

    Eigen::MatrixXd x1(dim,1);
    x1 << 0.1, 0.4;

    Eigen::MatrixXd x2(dim,1);
    x2 << 0.2, 0.7;
    
    Eigen::MatrixXd result = kernel2.Evaluate(x1,x2);

    Eigen::MatrixXd expected = A * kernel.Evaluate(x1,x2) * A.transpose();

    for(int j=0; j<A.rows(); ++j)
    {
	for(int i=0; i<A.rows(); ++i)
	    EXPECT_DOUBLE_EQ(expected(i,j), result(i,j));
    }
}



TEST(Approximation_GP, Clone)
{
    
    const unsigned dim = 2;
    auto kernel = ConstantKernel(dim, 2.0) * SquaredExpKernel(dim, 2.0, 0.35 );

    std::shared_ptr<KernelBase> kernel_copy = kernel.Clone();

    
    EXPECT_DOUBLE_EQ(kernel.inputDim, kernel_copy->inputDim);
    EXPECT_DOUBLE_EQ(kernel.coDim, kernel_copy->coDim);
    EXPECT_DOUBLE_EQ(kernel.numParams, kernel_copy->numParams);

    
    Eigen::VectorXd x1(dim);
    x1 << 0.1, 0.4;
    
    Eigen::VectorXd x2(dim);
    x2 << 0.2, 0.7;

    Eigen::MatrixXd result = kernel.Evaluate(x1,x2);
    Eigen::MatrixXd result_ptr = kernel_copy->Evaluate(x1,x2);
    
    EXPECT_DOUBLE_EQ(result(0,0), result_ptr(0,0));
}
