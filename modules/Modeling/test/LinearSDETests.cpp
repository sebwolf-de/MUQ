#include "gtest/gtest.h"

#include "MUQ/Modeling/LinearSDE.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Modeling;

TEST(LinearSDE, MeanCovariance)
{
    // Uses the Ornstein-Uhlenback SDE as a test case

    const int dim = 1;
    
    Eigen::MatrixXd F = -Eigen::MatrixXd::Identity(dim,dim);

    Eigen::MatrixXd L = Eigen::MatrixXd::Ones(dim, 1);
    
    Eigen::MatrixXd Q = Eigen::MatrixXd::Ones(1,1);

    boost::property_tree::ptree options;
    options.put("SDE.dt", 1e-3);

    LinearSDE sde(F,L,Q,options);
    
    Eigen::VectorXd mu0 = Eigen::VectorXd::Ones(dim);
    Eigen::MatrixXd cov0 = 0.5*Eigen::MatrixXd::Ones(dim,dim);

    Eigen::VectorXd mu;
    Eigen::MatrixXd cov;
    
    double endTime = 1.0;
    std::tie(mu, cov) = sde.EvolveDistribution(mu0, cov0, endTime);
    
    EXPECT_NEAR(mu0(0)*exp(F(0,0)*endTime), mu(0), 1e-3);
    EXPECT_NEAR(cov0(0,0), cov(0,0), 1e-3);
}
