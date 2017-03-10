#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;

TEST(Approximation_GP, HyperFit1)
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

    auto kernel = SquaredExpKernel(2.0,0.35, {0.1,10} )*PeriodicKernel(1.0,0.75,0.25, {0.5,5.0}, {0.5,5.0}, {0.25,0.5}) + WhiteNoiseKernel(1e-3, {1e-8,100});    
   
    // Create the GP
    ConstantMean mean(1);
    auto gp = ConstructGP(mean, kernel);


    gp.Fit(trainLocs, trainData);
    
    // Make a prediction
    auto post = gp.Predict(predLocs);

    //for(int i=0; i<predLocs.cols(); ++i)
    //	std::cout << "[" << predLocs(0,i) << ", " << post.mean(0,i) << "]," << std::endl;
    
}
