#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"
#include "MUQ/Utilities/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Utilities/LinearAlgebra/IdentityOperator.h"

#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"


#include "MUQ/Utilities/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Utilities/LinearAlgebra/LyapunovSolver.h"

#include "MUQ/Modeling/LinearSDE.h"


#include <unsupported/Eigen/Polynomials>
#include <Eigen/SparseCore>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace muq::Approximation;
using namespace muq::Utilities;


/** Implements the 2x2 block matrix in equation 28 of Solin and Sarkka */
class PeriodicKernel_F_block : public muq::Utilities::LinearOperator
{
public:
    PeriodicKernel_F_block(const double wjIn) : muq::Utilities::LinearOperator(2,2), wj(wjIn){};

    virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override
    {
        Eigen::MatrixXd output(2,x.cols());
        output.row(0) = -wj*x.row(1);
        output.row(1) = wj*x.row(0);
        return output;
    }
        
    virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override
    {
        Eigen::MatrixXd output(2,x.cols());
        output.row(0) = wj*x.row(1);
        output.row(1) = -wj*x.row(0);
        return output;
    }
    
    virtual Eigen::MatrixXd GetMatrix() override
    {
        Eigen::MatrixXd output(2,2);
        output << 0.0, -wj,
                   wj, 0.0;

        return output;
    }

private:
    const double wj;
    
};

std::shared_ptr<StateSpaceGP> PeriodicKernel::GetStateSpace(boost::property_tree::ptree sdeOptions) const
{

    // This is the same as J in the paper
    const int numTerms = sdeOptions.get("PeriodicKernel.StateSpace.NumTerms",4);
    
    const double w0 = 2.0*pi/period;

    const double l2 = std::pow(length, 2.0);

    
    Eigen::VectorXd q2s(numTerms+1); // holds all of the q_j^2 values from eqn 27 in "Explicit Link Between Periodic Covariance Functions and State Space Models"

    q2s(0) = boost::math::cyl_bessel_i(0, 1.0/l2) / exp(1.0/l2);
    for(int i=1; i<numTerms+1; ++i)
        q2s(i) = 2.0 * boost::math::cyl_bessel_i(i, 1.0/l2) / exp(1.0/l2);

    
    // SET UP F
    std::vector<std::shared_ptr<LinearOperator>> fBlocks(numTerms+1);
    for(int i=0; i<numTerms+1; ++i)
        fBlocks.at(i) = std::make_shared<PeriodicKernel_F_block>(w0*i);

    auto F = std::make_shared<BlockDiagonalOperator>(fBlocks);

    // SET UP L
    std::vector<std::shared_ptr<LinearOperator>> lBlocks(numTerms+1);
    for(int i=0; i<numTerms+1; ++i)
        lBlocks.at(i) = std::make_shared<IdentityOperator>(2);

    auto L = std::make_shared<BlockDiagonalOperator>(lBlocks);
    
    // Set up Pinf
    Eigen::MatrixXd Pinf = Eigen::MatrixXd::Zero(2*(numTerms+1), 2*(numTerms+1));
    for(int i=0; i<numTerms+1; ++i)
        Pinf.block(2*i, 2*i, 2, 2) = q2s(i) * Eigen::MatrixXd::Identity(2,2);

    // Set up Q
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(2*(numTerms+1), 2*(numTerms+1));
    
    // SET UP H
    std::vector<Eigen::Triplet<double>> Hcoeffs;
    for(int i=0; i<numTerms+1; ++i)
        Hcoeffs.push_back(Eigen::Triplet<double>(0, 2*i, 1.0));
        
    Eigen::SparseMatrix<double> Hmat(1, 2*(numTerms+1));
    Hmat.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());

    auto H = muq::Utilities::LinearOperator::Create(Hmat);

    // Create the SDE
    auto sde = std::make_shared<muq::Modeling::LinearSDE>(F,L,Q,sdeOptions);

    return std::make_shared<StateSpaceGP>(sde, H, Pinf);
    
}
