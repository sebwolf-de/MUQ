#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

#include <Eigen/Eigenvalues>

using namespace muq::Approximation;



KarhunenLoeveExpansion::KarhunenLoeveExpansion(std::shared_ptr<KernelBase> kernelIn,
                                                  Eigen::MatrixXd      const& seedPtsIn,
                                                  Eigen::VectorXd      const& seedWtsIn,
                                                  boost::property_tree::ptree options) : seedPts(seedPtsIn), seedWts(seedWtsIn), covKernel(kernelIn)
{
    // Get the type of truncation to use (Energy, FixedNumber)
    std::string truncType = options.get("KarhunenLoeve.TruncationType","Energy");
    if( truncType.compare("Energy") && truncType.compare("FixedNumber") )
    {
      std::cerr << "\nERROR:  Invalid options for \"KarhunenLoeve.TruncationType\" set in KarhunenLoveExpansion::KarhunenLoeveExpansion.  Valide options are \"Energy\" or \"FixedNumber\", but a value of " + truncType + " was given.\n\n";
      assert(!truncType.compare("Energy") || !truncType.compare("FixedNumber"));
    }

    // We will approximation the KL modes as the GP with known values at the seed points.  To get those values, we first need to solve the discrete eigenvalue problem
    Eigen::MatrixXd seedCov = covKernel->BuildCovariance(seedPts);
    
    // Scale the columns of seedCov by the seed weights
    for(int i=0; i<seedWts.size(); ++i)
      seedCov.col(i) *= seedWts(i);

    Eigen::EigenSolver<Eigen::MatrixXd> eigSolver;
    eigSolver.compute(seedCov);

    // Are the eigenvalues increasing or decreasing...
    bool eigsIncrease = eigSolver.eigenvalues()(0).real() < eigSolver.eigenvalues()(seedWts.size()-1).real();

    // Figure out how many modes we want to keep
    int numModes;
    if(!truncType.compare("FixedNumber")){
      numModes = options.get("KarhunenLoeve.NumModes", seedWts.size());

    }else if(!truncType.compare("Energy")){

      double energyFraction = options.get("KarhunenLoeve.EnergyTol", 0.9);
      assert(energyFraction>0);
      assert(energyFraction<=1);

      double totalEnergy = eigSolver.eigenvalues().squaredNorm();
      numModes = 0;
      int ind;
      double cumEnergy=0;
      while((cumEnergy < energyFraction * totalEnergy)&&(numModes<seedWts.size()))
      {
        if(eigsIncrease){
          ind = seedWts.size()-1-numModes;
        }else{
          ind = numModes;
        }

        cumEnergy += std::pow(eigSolver.eigenvalues()(ind).real(), 2.0);
        ++numModes;
      }
    }

    if(eigsIncrease){
      modeEigs = eigSolver.eigenvalues().tail(numModes).real().reverse();
      modeVecs = eigSolver.eigenvectors().rightCols(numModes).real().rowwise().reverse();
    }else{
      modeEigs = eigSolver.eigenvalues().tail(numModes).real();
      modeVecs = eigSolver.eigenvectors().rightCols(numModes).real();
    }

    modeScales = 1.0/modeEigs.array();
}

Eigen::MatrixXd KarhunenLoeveExpansion::GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts)
{

    // Build the cross covariance between the seed points and the evaluation points
    Eigen::MatrixXd crossCov = covKernel->BuildCovariance(pts,seedPts);

    
    return crossCov * seedWts.asDiagonal() * modeVecs * modeScales.asDiagonal();
}

Eigen::VectorXd KarhunenLoeveExpansion::Evaluate(Eigen::Ref<const Eigen::MatrixXd> const& pts,
                                                 Eigen::Ref<const Eigen::VectorXd> const& coeffs)
{
    return GetModes(pts) * (modeEigs.array().sqrt()*coeffs.array()).matrix();
}
