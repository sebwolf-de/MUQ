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


TEST(LinearSDE, Concatenate)
{
    // Uses the Ornstein-Uhlenback SDE as a test case
    const int dim = 1;
    
    Eigen::MatrixXd F = -Eigen::MatrixXd::Identity(dim,dim);

    Eigen::MatrixXd L = Eigen::MatrixXd::Ones(dim, 1);
    
    Eigen::MatrixXd Q = Eigen::MatrixXd::Ones(1,1);

    boost::property_tree::ptree options;
    options.put("SDE.dt", 1e-3);

    std::vector<std::shared_ptr<LinearSDE>> sdes(2);
    sdes.at(0) = std::make_shared<LinearSDE>(F,L,Q,options);
    sdes.at(1) = std::make_shared<LinearSDE>(F,L,Q,options);

    auto sde = LinearSDE::Concatenate(sdes, options);

    EXPECT_EQ(2*dim, sde->stateDim);
    EXPECT_EQ(2*dim, sde->GetF()->rows());
    EXPECT_EQ(2*dim, sde->GetF()->cols());
    EXPECT_EQ(2*dim, sde->GetL()->rows());
    EXPECT_EQ(2*dim, sde->GetL()->rows());
    
    Eigen::VectorXd mu0 = Eigen::VectorXd::Ones(2*dim);
    Eigen::MatrixXd cov0 = 0.5*Eigen::MatrixXd::Identity(2*dim,2*dim);

    Eigen::VectorXd mu;
    Eigen::MatrixXd cov;
    
    double endTime = 1.0;
    std::tie(mu, cov) = sde->EvolveDistribution(mu0, cov0, endTime);

    EXPECT_NEAR(mu0(0)*exp(F(0,0)*endTime), mu(0), 1e-3);
    EXPECT_NEAR(mu0(1)*exp(F(0,0)*endTime), mu(1), 1e-3);
    EXPECT_NEAR(cov0(0,0), cov(0,0), 1e-3);
    EXPECT_NEAR(cov0(0,1), cov(0,1), 1e-3);
    EXPECT_NEAR(cov0(1,1), cov(1,1), 1e-3);
    EXPECT_NEAR(cov0(1,0), cov(1,0), 1e-3);
}
