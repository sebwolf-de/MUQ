#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"
#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Utilities/Exceptions.h"

#include <Eigen/Dense>

using namespace muq::Approximation;
using namespace muq::Modeling;
using namespace muq::Utilities;

StateSpaceGP::StateSpaceGP(std::shared_ptr<MeanFunctionBase> meanIn,
                           std::shared_ptr<KernelBase>       covKernelIn,
                           boost::property_tree::ptree       options) : StateSpaceGP(covKernelIn->GetStateSpace(options),
                                                                                     meanIn,
                                                                                     covKernelIn)
{
};
    

StateSpaceGP::StateSpaceGP(std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> ssInfo,
                           std::shared_ptr<MeanFunctionBase> meanIn,
                           std::shared_ptr<KernelBase>       covKernelIn) : GaussianProcess(meanIn, covKernelIn), stateDim(std::get<2>(ssInfo).rows()),
                                                              sde(std::get<0>(ssInfo)),
                                                              obsOp(std::get<1>(ssInfo)),
                                                              L(std::get<2>(ssInfo).selfadjointView<Eigen::Lower>().llt().matrixL()),
                                                              mean(meanIn),
                                                              covKernel(covKernelIn)
{
}

Eigen::MatrixXd StateSpaceGP::Sample(Eigen::MatrixXd const& times)
{

    // Generate sample for initial condition
    Eigen::VectorXd x = L.triangularView<Eigen::Lower>()*RandomGenerator::GetNormal(L.rows());
    
    // Make space for the simulated GP
    Eigen::MatrixXd output(obsOp->rows(), times.size());

    output.col(0) = obsOp->Apply(x);

    // Step through the each time and integrate the SDE between times
    for(int i=0; i<times.size()-1; ++i)
    {
        x = sde->EvolveState(x, times(i+1)-times(i));
        output.col(i+1) = obsOp->Apply(x);
    }

    return output;
}


Eigen::MatrixXd StateSpaceGP::PredictMean(Eigen::MatrixXd const& newPts)
{

    return Eigen::MatrixXd::Zero(coDim, newPts.size());
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> StateSpaceGP::Predict(Eigen::MatrixXd const& times,
                                                                  CovarianceType         covType)
{

    std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::MatrixXd>> output;
    output.first.resize(times.cols());
    output.second.resize(times.cols());
    
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> currDist;
    currDist.second  = L.triangularView<Eigen::Lower>()*L.transpose();
    currDist.first = Eigen::VectorXd::Zero(stateDim);

    output.first.at(0)  = obsOp->Apply(currDist.first);
    output.second.at(0) = obsOp->Apply(obsOp->Apply(currDist.second).transpose());

    // Step through the each time and integrate the SDE between times
    for(int i=1; i<times.cols(); ++i)
    {
        currDist = sde->EvolveDistribution(currDist, times(i)-times(i-1));
        
        output.first.at(i) = obsOp->Apply(currDist.first);
        output.second.at(i) = obsOp->Apply(obsOp->Apply(currDist.second).transpose());
    }

    return output;
}


double StateSpaceGP::LogLikelihood(Eigen::MatrixXd const& xs,
                                   Eigen::MatrixXd const& vals)
{
    return 0.0;
}

double StateSpaceGP::MarginalLogLikelihood(Eigen::Ref<Eigen::VectorXd> grad,
                                           bool                        computeGrad)
{
    return 0.0;
}





void StateSpaceGP::SetObs(std::shared_ptr<muq::Utilities::LinearOperator> newObs)
{

    if(newObs->cols() != obsOp->cols())
        throw muq::WrongSizeError("In StateSpaceGP::SetObs: The new observation operator has " + std::to_string(newObs->cols()) + " columns, which does not match the system dimension " + std::to_string(obsOp->cols()));

    obsOp = newObs;
}

std::shared_ptr<StateSpaceGP> StateSpaceGP::Concatenate(std::vector<std::shared_ptr<StateSpaceGP>> const& gps)
{
    // Build the concatenated SDE
    std::vector<std::shared_ptr<muq::Modeling::LinearSDE>> sdes(gps.size());
    for(int i=0; i<gps.size(); ++i)
        sdes.at(i) = gps.at(i)->GetSDE();

    auto sde = LinearSDE::Concatenate(sdes);
    
    // Build a concatenated observation operator
    std::vector<std::shared_ptr<muq::Utilities::LinearOperator>> obsOps(gps.size());
    for(int i=0; i<gps.size(); ++i)
        obsOps.at(i) = gps.at(i)->GetObs();
    
    auto H = std::make_shared<BlockDiagonalOperator>(obsOps);

    // Build the concatenated covariance
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(sde->stateDim, sde->stateDim);
    int currRow = 0;
    for(int i=0; i<gps.size(); ++i)
    {
        Q.block(currRow,currRow, gps.at(i)->stateDim, gps.at(i)->stateDim) = gps.at(i)->GetCov();
        currRow += gps.at(i)->stateDim;
    }
  
    return std::make_shared<StateSpaceGP>(sde, H, Q);
}

                 
