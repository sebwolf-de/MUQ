#ifndef STATESPACEGP_H
#define STATESPACEGP_H

#include <memory>

#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearSDE.h"

#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include <Eigen/Core>

namespace muq
{
namespace Approximation
{


class StateSpaceGP : public GaussianProcess
{
public:

    StateSpaceGP(MeanFunctionBase&           meanIn,
                 KernelBase&                 kernelIn,
                 boost::property_tree::ptree options = boost::property_tree::ptree()) : StateSpaceGP(meanIn.Clone(), kernelIn.Clone(), options){};
    
    StateSpaceGP(std::shared_ptr<MeanFunctionBase> meanIn,
                 std::shared_ptr<KernelBase>       covKernelIn,
                 boost::property_tree::ptree       options = boost::property_tree::ptree());

        
    virtual Eigen::MatrixXd Sample(Eigen::MatrixXd const& times) override;
    
    
    virtual std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Predict(Eigen::MatrixXd const& newLocs,
                                                                CovarianceType         covType) override;

    virtual Eigen::MatrixXd PredictMean(Eigen::MatrixXd const& newPts) override;

    
    virtual double LogLikelihood(Eigen::MatrixXd const& xs,
                                 Eigen::MatrixXd const& vals) override;
    
    virtual double MarginalLogLikelihood(Eigen::Ref<Eigen::VectorXd> grad,
                                         bool                        computeGrad = true) override;

    
    std::shared_ptr<muq::Modeling::LinearSDE> GetSDE(){return sde;};
    
    std::shared_ptr<muq::Utilities::LinearOperator> GetObs(){return obsOp;};

    void SetObs(std::shared_ptr<muq::Utilities::LinearOperator> newObs);

    Eigen::MatrixXd GetCov(){return L.triangularView<Eigen::Lower>()*L.transpose();};


    
    const int stateDim;

    template<typename KernelType1, typename...  KernelTypes>
    static std::shared_ptr<StateSpaceGP> Concatenate(KernelType1 const& k1, KernelTypes... kernels)
    {
      std::vector<std::shared_ptr<KernelBase>> kernelVec(1);
      kernelVec.at(0) = k1.Clone();
      return Concatenate(kernelVec, kernels...);
    }

    template<typename KernelType1, typename...  KernelTypes>
    static std::shared_ptr<StateSpaceGP> Concatenate(std::vector<std::shared_ptr<KernelBase>> &kernelVec, 
                                                     KernelType1 const& k1,
                                                     KernelTypes... kernels)
    {
      kernelVec.push_back(k1.Clone());
      return Concatenate(kernelVec, kernels...);
    }

    template<typename KernelType1>
    static std::shared_ptr<StateSpaceGP> Concatenate(std::vector<std::shared_ptr<KernelBase>> &kernelVec,
                                                       KernelType1 const& k1)
    {
        kernelVec.push_back(k1);
        return Concatenate(kernelVec);
    }

    static std::shared_ptr<StateSpaceGP> Concatenate(std::vector<std::shared_ptr<KernelBase>> const& gps);
    static std::shared_ptr<StateSpaceGP> Concatenate(std::vector<std::shared_ptr<StateSpaceGP>> const& gps);
    
protected:

    StateSpaceGP(std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> ssInfo,
                 std::shared_ptr<MeanFunctionBase> meanIn,
                 std::shared_ptr<KernelBase>       covKernelIn);


    void SortObservations();

    bool ComputeAQ(double dt);
    Eigen::MatrixXd sdeA, sdeQ;
    double dtAQ; // the last deltat passed to the ComputeAQ function
    
    // Stocastic Differential equation describing correlations
    std::shared_ptr<muq::Modeling::LinearSDE> sde;

    // Observation operator
    std::shared_ptr<muq::Utilities::LinearOperator> obsOp;
    
    // Cholesky factor of the stationary covariance matrix (used to initialize SDE integration)
    Eigen::MatrixXd L;

    std::shared_ptr<MeanFunctionBase> mean;
    std::shared_ptr<KernelBase>       covKernel;
    

};


} // namespace Approximation
} // namespace GP





#endif
