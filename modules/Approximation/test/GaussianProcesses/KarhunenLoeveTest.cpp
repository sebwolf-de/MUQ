
#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/Exceptions.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Utilities;
using namespace muq::Approximation;

TEST(Approximation_GP, KarhunenLoeve_BasicConstruction)
{

    double sigma2 = 1.0;
    double length = 0.2;
    double nu = 3.0/2.0;

    auto kernel = std::make_shared<MaternKernel>(1, sigma2, length, nu);

    const int numSeeds = 20;
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
    std::cout << "Number of modes = " << modes.cols() << std::endl;

    const unsigned int N = 1e4;

    // Picke a row of the modes and evaluate the variance at that point
    Eigen::VectorXd samps = modes.row(1) * RandomGenerator::GetNormal(modes.cols(), N);
    double mu = samps.mean();
    EXPECT_NEAR(mu, 0.0, 1.0e-3);

    double var = (samps-Eigen::VectorXd::Constant(N,mu)).array().pow(2.0).sum()/(N-1.0);
    EXPECT_NEAR(var, sigma2, 1.0e-3);

    // check the mean
    const double x = 0.25;
    for( unsigned int i=0; i<N; ++i )
      samps(i) = kl.Evaluate(Eigen::VectorXd::Constant(1, x), RandomGenerator::GetNormal(modes.cols()))(0);

    mu = samps.mean();
    EXPECT_NEAR(mu, 0.0, 1.0e-3);

    // check the variance
    var = (samps-Eigen::VectorXd::Constant(N,mu)).array().pow(2.0).sum()/(N-1.0);

    EXPECT_NEAR(var, sigma2, 1.0e-3);
}
