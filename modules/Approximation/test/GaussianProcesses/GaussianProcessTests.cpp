#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;
using namespace muq::Utilities;

TEST(Approximation_GP, HyperFit1d)
{
    const unsigned numPred  = 50;
    const unsigned maxTrain = 10;
    
    const double pi = 4.0 * atan(1.0);

    // Set linearly space locations where we want to evaluate the Gaussian Process
    Eigen::RowVectorXd predLocs(numPred);
    predLocs.setLinSpaced(numPred, 0,1);

    // Generate random training locations
    Eigen::RowVectorXd trainLocs = Eigen::RowVectorXd::LinSpaced(maxTrain, 0, 1);
    
    // Generate the training data (with some random noise)
    Eigen::RowVectorXd trainData(maxTrain);
    for(int i=0; i<maxTrain; ++i)
	trainData(i) = sin(4*2*pi*trainLocs(i) );

    trainData += sqrt(1e-4)*RandomGenerator::GetNormal(maxTrain).transpose();
    
    const unsigned dim = 1;
    auto kernel = SquaredExpKernel(dim, 2.0, 0.35, {0.1,10} ) * PeriodicKernel(dim, 1.0, 0.75, 0.25, {0.5,5.0}, {0.5,5.0}, {0.25,0.5}) + WhiteNoiseKernel(dim, 1e-3, {1e-8,100});    
   
    // Create the GP
    ConstantMean mean(dim, 1);
    GaussianProcess gp(mean, kernel);


    for(int i=0; i<trainData.size(); ++i)
        gp.Condition(trainLocs.col(i), trainData.col(i));

    // Fit the hyperparameters
    gp.Optimize();
    
    // Make a prediction
    Eigen::MatrixXd postMean = gp.PredictMean(predLocs);

    //for(int i=0; i<predLocs.cols(); ++i)
    //	std::cout << "[" << predLocs(0,i) << ", " << postMean(0,i) << "]," << std::endl;
    
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

    for(int i=0; i<trainLocs.cols(); ++i)
        gp.Condition(trainLocs.col(i),trainData.col(i));

    gp.Optimize();
    
    // Make a prediction
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> post = gp.Predict(predLocs, GaussianProcess::FullCov);

}


