#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include "MUQ/Utilities/Exceptions.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;






TEST(Approximation_GP, MaternStateSpace)
{

    const double sigma2 = 1.0;
    const double length = 0.15;

    const double nu = 3.0/2.0;
    
    MaternKernel kernel(1, sigma2, length, nu);

    ConstantMean mu(1,1);
    StateSpaceGP gp(mu, kernel);
    
    EXPECT_EQ(nu+0.5, gp.stateDim);

    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(100, 0, 1);
    Eigen::MatrixXd realization = gp.Sample(obsTimes);

}

TEST(Approximation_GP, PeriodicStateSpace)
{

    const double sigma2 = 1.0;
    const double length = 0.6;
    const double period = 0.25;
    const double periodN = 50; // how many steps per period
    
    PeriodicKernel kernel(1, sigma2, length, period);

    boost::property_tree::ptree options;
    options.put("PeriodicKernel.StateSpace.NumTerms",8);
    options.put("SDE.dt", 5e-5);

    ConstantMean mu(1,1);
    StateSpaceGP gp(mu, kernel, options);
    
    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(5*periodN+1, 0, 5*period);
        
    Eigen::MatrixXd realization = gp.Sample(obsTimes);

    // Make sure the sample is periodic
    for(int i=0; i<obsTimes.size()-periodN-1; ++i)
        EXPECT_NEAR(realization(i), realization(i+periodN), 1e-1);
    
}


TEST(Approximation_GP, ProductStateSpace)
{

    const double sigma2 = 1.0;
    const double length = 0.6;
    const double nu = 3.0/2.0;
    const double period = 0.25;
    const double periodN = 50; // how many steps per period
    
    PeriodicKernel kernel1(1, sigma2, 0.8, period);
    MaternKernel kernel2(1, sigma2, 2.0, nu);

    auto kernel12 = kernel1*kernel2;
    auto kernel21 = kernel2*kernel1;
    auto kernel22 = kernel2*kernel2;

    boost::property_tree::ptree options;
    options.put("PeriodicKernel.StateSpace.NumTerms",7);
    options.put("SDE.dt", 1e-4);

    ConstantMean mu(1,1);
    StateSpaceGP gp1(mu, kernel1, options);
    StateSpaceGP gp2(mu, kernel2, options);
    EXPECT_EQ(int(nu+0.5), gp2.stateDim);

    StateSpaceGP gp12(mu, kernel12, options);
    StateSpaceGP gp22(mu, kernel21, options);
    
    EXPECT_THROW(kernel22.GetStateSpace(options), muq::NotImplementedError);
    EXPECT_EQ(gp1.stateDim *gp2.stateDim, gp12.stateDim);
    
    
    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(5*periodN+1, 0, 5*period);

    Eigen::MatrixXd realization1  = gp1.Sample(obsTimes);
    Eigen::MatrixXd realization2  = gp2.Sample(obsTimes);
    Eigen::MatrixXd realization12 = gp12.Sample(obsTimes);
    
}
