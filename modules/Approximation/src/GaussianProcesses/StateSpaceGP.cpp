#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <Eigen/Dense>

using namespace muq::Approximation;
using namespace muq::Modeling;
using namespace muq::Utilities;


StateSpaceGP::StateSpaceGP(std::shared_ptr<LinearSDE>      sdeIn,
                           std::shared_ptr<LinearOperator> obsOpIn,
                           Eigen::MatrixXd          const& pInfIn) : stateDim(pInfIn.rows()), sde(sdeIn), obsOp(obsOpIn), L(pInfIn.selfadjointView<Eigen::Lower>().llt().matrixL())
{

}


Eigen::VectorXd StateSpaceGP::Sample(Eigen::VectorXd const& times) const
{

    // Generate sample for initial condition
    Eigen::VectorXd x = L.triangularView<Eigen::Lower>()*RandomGenerator::GetNormal(L.rows());
    std::cout << "x = " << x << std::endl;
    
    // Make space for the simulated GP
    Eigen::VectorXd output(times.size());

    output(0) = obsOp->Apply(x)(0);
    
    // Step through the each time and integrate the SDE between times
    for(int i=0; i<times.size()-1; ++i)
    {   
        x = sde->EvolveState(x, times(i+1)-times(i));
        output(i+1) = obsOp->Apply(x)(0);
    }

    return output;
}


std::pair<Eigen::VectorXd, Eigen::VectorXd> StateSpaceGP::ComputeMeanVar(Eigen::VectorXd const& times) const
{

    Eigen::VectorXd outputMean(times.size());
    Eigen::VectorXd outputVar(times.size());

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> currDist;
    currDist.second  = L.triangularView<Eigen::Lower>()*L.transpose();
    currDist.first = Eigen::VectorXd::Zero(stateDim);

    outputMean(0) = obsOp->Apply(currDist.first)(0);
    outputVar(0)  = obsOp->Apply(obsOp->Apply(currDist.second).transpose())(0,0);
    
    // Step through the each time and integrate the SDE between times
    for(int i=1; i<times.size(); ++i)
    {
        currDist = sde->EvolveDistribution(currDist, times(i+1)-times(i));

        outputMean(i) = obsOp->Apply(currDist.first)(0);
        outputVar(i) = obsOp->Apply(obsOp->Apply(currDist.second).transpose())(0,0);
    }

    return std::make_pair(outputMean, outputVar);
}




                 
