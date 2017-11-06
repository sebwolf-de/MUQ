#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

#include "MUQ/Utilities/Exceptions.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;


TEST(Approximation_GP, KarhunenLoeve_BasicConstruction)
{

    double sigma2 = 1.0;
    double length = 0.2;
    double nu = 3.0/2.0;
    
    auto kernel = std::make_shared<MaternKernel>(1, sigma2, length, nu);

    const int numSeeds = 10;
    const int numEval = 200;
    
    const double lb = 0.0;
    const double ub = 1.0;
    
    Eigen::MatrixXd seedPts(1,numSeeds);
    for(int i=0; i<numSeeds; ++i)
        seedPts(0, i) = lb + (ub-lb)*double(i)/double(numSeeds);

    Eigen::VectorXd seedWts = (ub-lb)/double(numSeeds)*Eigen::VectorXd::Ones(numSeeds);

    // Construct the Karhunen-Loeve decomposition based on a discrete eigenvalue problem at the seed points
    KarhunenLoeveExpansion kl(kernel, seedPts, seedWts);

    // Evaluate the KL expansion at a bunch of other points
    Eigen::MatrixXd evalPts(1,numEval);
    for(int i=0; i<numEval; ++i)
        evalPts(0, i) = lb + (ub-lb)*double(i)/double(numEval);
    
    Eigen::MatrixXd modes = kl.GetModes(evalPts);
   
}
