#ifndef STATESPACEGP
#define STATESPACEGP

#include <memory>

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearSDE.h"

#include <Eigen/Core>

namespace muq
{
namespace Approximation
{

    
class StateSpaceGP
{
public:
    StateSpaceGP(std::shared_ptr<muq::Modeling::LinearSDE>       sdeIn,
                 std::shared_ptr<muq::Utilities::LinearOperator> obsOpIn,
                 Eigen::MatrixXd                          const& pInfIn);
                 

    Eigen::VectorXd Sample(Eigen::VectorXd const& times) const;

    std::pair<Eigen::VectorXd, Eigen::VectorXd> ComputeMeanVar(Eigen::VectorXd const& times) const;

    
    std::shared_ptr<muq::Modeling::LinearSDE> GetSDE(){return sde;};
    
    std::shared_ptr<muq::Utilities::LinearOperator> GetObs(){return obsOp;};

    Eigen::MatrixXd GetCov(){return L.triangularView<Eigen::Lower>()*L.transpose();};

    const int stateDim;

    static std::shared_ptr<StateSpaceGP> Concatenate(std::vector<std::shared_ptr<StateSpaceGP>> const& gps);
    
private:

    // Stocastic Differential equation describing correlations
    std::shared_ptr<muq::Modeling::LinearSDE> sde;

    // Observation operator
    std::shared_ptr<muq::Utilities::LinearOperator> obsOp;
    
    // Cholesky factor of the stationary covariance matrix (used to initialize SDE integration)
    Eigen::MatrixXd L; 

};


} // namespace Approximation
} // namespace GP





#endif
