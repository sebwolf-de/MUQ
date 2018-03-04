#ifndef SUMKERNEL_H
#define SUMKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


namespace muq
{
namespace Approximation
{

/**

@class SumKernel
@ingroup CovarianceKernels
\f[
k)x,y) = k_1(x,y) + k_2(x,y)
\f]

 */
class SumKernel : public KernelBase
{

public:
    SumKernel(std::shared_ptr<KernelBase> kernel1In,
              std::shared_ptr<KernelBase> kernel2In) : KernelBase(kernel1In->inputDim,
                                    														  kernel1In->coDim,
                                    														  kernel1In->numParams + kernel2In->numParams),
	                                                     kernel1(kernel1In),
	                                                     kernel2(kernel2In)
    {
      assert(kernel2->inputDim == kernel2->inputDim);
	    assert(kernel1->coDim == kernel2->coDim);
    };

    virtual ~SumKernel(){};


    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override
    {

        kernel1->FillBlock(x1,x2, params.head(kernel1->numParams), block);

	      Eigen::MatrixXd temp(coDim, coDim);
        kernel2->FillBlock(x1,x2, params.tail(kernel2->numParams), temp);

        block += temp;
    }


    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const override
    {
        std::vector<std::shared_ptr<muq::Modeling::LinearSDE>> sdes(2);
        std::vector<std::shared_ptr<muq::Utilities::LinearOperator>> linOps(2);
        Eigen::MatrixXd pinf1, pinf2;

        // Get the statespace information from each component
        std::tie(sdes.at(0), linOps.at(0), pinf1) = kernel1->GetStateSpace(sdeOptions);
        std::tie(sdes.at(1), linOps.at(1), pinf2) = kernel2->GetStateSpace(sdeOptions);

        // Concantenate the sdes
        auto newSDE = muq::Modeling::LinearSDE::Concatenate(sdes,sdeOptions);

        // Concatenate the linear operators
        auto newObsOp = std::make_shared<muq::Utilities::BlockRowOperator>(linOps);

        // Set up the combined stationary covariance
        Eigen::MatrixXd newPinf = Eigen::MatrixXd::Zero(pinf1.rows() + pinf2.rows(), pinf1.cols() + pinf2.cols());
        newPinf.block(0,0,pinf1.rows(),pinf1.cols()) = pinf1;
        newPinf.block(pinf1.rows(),pinf1.cols(), pinf2.rows(), pinf2.cols()) = pinf2;

        return std::make_tuple(newSDE, newObsOp, newPinf);

    }

private:
    std::shared_ptr<KernelBase>  kernel1;
    std::shared_ptr<KernelBase> kernel2;
};

}
}


#endif
