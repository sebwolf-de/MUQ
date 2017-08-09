#include "MUQ/Modeling/LinearSDE.h"

#include <random>

#include <unsupported/Eigen/MatrixFunctions>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"

using namespace muq::Modeling;
using namespace muq::Utilities;

LinearSDE::LinearSDE(std::shared_ptr<LinearOperator>    Fin,
                     std::shared_ptr<LinearOperator>    Lin,
                     Eigen::MatrixXd             const& Qin,
                     boost::property_tree::ptree        options) : stateDim(Fin->rows()), F(Fin), L(Lin), Q(Qin)
{
    // Extract options from the ptree
    ExtractOptions(options);
    
    // Compute the Cholesky decomposition of the white noise process
    sqrtQ = Q.llt().matrixL();
};

void LinearSDE::ExtractOptions(boost::property_tree::ptree options)
{
    dt = options.get("SDE.dt",1e-4);
}

Eigen::VectorXd LinearSDE::EvolveState(Eigen::VectorXd const& f0,
                                       double                 T) const
{
    Eigen::VectorXd f = f0;

    const int numTimes = std::ceil(T/dt);

    Eigen::VectorXd z;
    
    // Take all but the last step.  The last step might be a partial step
    for(int i=0; i<numTimes-1; ++i)
    {
        z = sqrt(dt) * (sqrtQ.triangularView<Eigen::Lower>() * RandomGenerator::GetNormal(sqrtQ.cols()) ).eval();
        
        f += dt*F->Apply(f) + L->Apply( z );
    }

    // Now take the last step
    double lastDt = T-(numTimes-1)*dt;
    
    z = sqrt(lastDt) * (sqrtQ.triangularView<Eigen::Lower>() * RandomGenerator::GetNormal(sqrtQ.cols())).eval();

    f += lastDt*F->Apply(f) + L->Apply( z );
    
    return f;
}


std::pair<Eigen::VectorXd, Eigen::MatrixXd> LinearSDE::EvolveDistribution(Eigen::VectorXd const& mu0,
                                                                          Eigen::MatrixXd const& gamma0,
                                                                          double                 T) const
{

    Eigen::VectorXd mu = mu0;
    Eigen::MatrixXd gamma = gamma0;

    const int numTimes = std::ceil(T/dt);

    Eigen::MatrixXd LQLT = L->Apply( L->Apply(Q).transpose().eval() );
    LQLT = 0.5*(LQLT + LQLT.transpose()); // <- Make sure LQLT is symmetric
    
    Eigen::MatrixXd Fgamma;

    // Take all but the last step because the last step might be a partial step.
    for(int i=0; i<numTimes-1; ++i)
    {
        mu += dt * F->Apply(mu);
        Fgamma = F->Apply(gamma);
        gamma += dt * (Fgamma + Fgamma.transpose() + LQLT);
    }

    // Take the last step
    double lastDt = T-(numTimes-1)*dt;
    mu += lastDt * F->Apply(mu);
    Fgamma = F->Apply(gamma);
    gamma += lastDt * (Fgamma + Fgamma.transpose() + LQLT);

    return std::make_pair(mu,gamma);
}


std::shared_ptr<LinearSDE> LinearSDE::Concatenate(std::vector<std::shared_ptr<LinearSDE>> const& sdes,
                                                  boost::property_tree::ptree                    options)
{

    int stateDim = 0;
    int stochDim = 0;
    for(auto& sde : sdes){
        stateDim += sde->stateDim;
        stochDim += sde->L->cols();
    }
    
    std::vector<std::shared_ptr<LinearOperator>> Fs(sdes.size());
    std::vector<std::shared_ptr<LinearOperator>> Ls(sdes.size());

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(stochDim, stochDim);

    int currDim = 0;
    for(int i=0; i<sdes.size(); ++i){
        Fs.at(i) = sdes.at(i)->GetF();
        Ls.at(i) = sdes.at(i)->GetL();

        Q.block(currDim, currDim, Ls.at(i)->cols(), Ls.at(i)->cols()) = sdes.at(i)->GetQ();

        currDim += Ls.at(i)->cols();
    }

    auto F = std::make_shared<BlockDiagonalOperator>(Fs);
    auto L = std::make_shared<BlockDiagonalOperator>(Ls);

    return std::make_shared<LinearSDE>(F,L,Q, options);
}
