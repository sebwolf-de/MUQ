#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

#include <Eigen/Eigenvalues>

using namespace muq::Approximation;



KarhunenLoeveExpansion::KarhunenLoeveExpansion(std::shared_ptr<KernelBase> kernelIn,
                                                  Eigen::MatrixXd      const& seedPtsIn,
                                                  Eigen::VectorXd      const& seedWtsIn,
                                                  boost::property_tree::ptree options) : seedPts(seedPtsIn), seedWts(seedWtsIn), covKernel(kernelIn)
{

    int numModes = options.get("KarhunenLoeve.NumModes", seedWts.size());
    
    // We will approximation the KL modes as the GP with known values at the seed points.  To get those values, we first need to solve the discrete eigenvalue problem
    Eigen::MatrixXd seedCov = covKernel->BuildCovariance(seedPts);
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver;
    eigSolver.compute(seedCov);

    double minEig = eigSolver.eigenvalues().minCoeff();
    if(minEig<=0){
        std::cerr << "YIKES!  The covariance matrix wasn't positive definite.";
        assert(minEig>0);
    }
    
    modeEigs = eigSolver.eigenvalues().tail(numModes).reverse();
    modeVecs = eigSolver.eigenvectors().rightCols(numModes).rowwise().reverse();
}

Eigen::MatrixXd KarhunenLoeveExpansion::GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts)
{

    // Build the cross covariance between the seed points and the evaluation points
    Eigen::MatrixXd crossCov = covKernel->BuildCovariance(pts,seedPts);

    Eigen::VectorXd scale = seedWts.array()/modeEigs.array();
    return crossCov * modeVecs * scale.asDiagonal();
}

Eigen::VectorXd KarhunenLoeveExpansion::Evaluate(Eigen::Ref<const Eigen::MatrixXd> const& pts,
                                                 Eigen::Ref<const Eigen::VectorXd> const& coeffs)
{
    return GetModes(pts) * (modeEigs.array().sqrt()*coeffs.array()).matrix();
}
