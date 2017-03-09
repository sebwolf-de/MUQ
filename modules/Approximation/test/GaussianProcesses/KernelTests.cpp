#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>


TEST(Approximation_GP, HyperFit1)
{
    const unsigned numPred  = 200;
    const unsigned maxTrain = 1000;
    
    const double pi = 4.0 * atan(1.0); //boost::math::constants::pi<double>();
    
    Eigen::RowVectorXd predLocs(numPred);
    predLocs.setLinSpaced(numPred, 0,1);

    Eigen::RowVectorXd trainLocs(maxTrain);
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

}
