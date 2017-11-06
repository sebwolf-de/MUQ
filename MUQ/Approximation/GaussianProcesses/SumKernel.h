#ifndef SUMKERNEL_H
#define SUMKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


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

template<typename LeftType, typename RightType>
class SumKernel : public KernelImpl<SumKernel<LeftType,RightType>>
{

public:
    SumKernel(LeftType const& kernel1In, RightType const& kernel2In) : KernelImpl<SumKernel<LeftType,RightType>>(kernel1In.inputDim,
														 kernel1In.coDim,
														 kernel1In.numParams + kernel2In.numParams),
	                                                               kernel1(kernel1In),
	                                                               kernel2(kernel2In)
    {
        assert(kernel2.inputDim == kernel2.inputDim);
	assert(kernel1.coDim == kernel2.coDim);
    };

    virtual ~SumKernel(){};
    
    template<typename VecType, typename MatType>
    inline void EvaluateImpl(VecType const& x1, VecType  const& x2 , MatType & cov) const
    {
	Eigen::MatrixXd temp1(this->coDim, this->coDim);
	Eigen::MatrixXd temp2(this->coDim, this->coDim);
	
        kernel1.EvaluateImpl(x1,x2, temp1);
	kernel2.EvaluateImpl(x1,x2, temp2);
	cov = temp1 + temp2;
    }

    template<typename VecType1, typename VecType2, typename MatType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < this->numParams);

        if(wrt < kernel1.numParams )
	{
	    return kernel1.GetDerivative(x1, x2, wrt, derivs);
	}
	else
	{
	    return kernel2.GetDerivative(x1,x2, wrt-kernel1.numParams, derivs);
	}
    }
    
    virtual Eigen::MatrixXd GetParamBounds() const override
    {
	Eigen::MatrixXd bounds(2,this->numParams);

	bounds.block(0,0,2,kernel1.numParams) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.numParams,2,kernel2.numParams) = kernel2.GetParamBounds();

	return bounds;
    };
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(this->numParams);

	params.head(kernel1.numParams) = kernel1.GetParams();
	params.tail(kernel2.numParams) = kernel2.GetParams();

	return params;
    };
        
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        kernel1.SetParams(params.head(kernel1.numParams));
	kernel2.SetParams(params.tail(kernel2.numParams));
    }

    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const override
    {
        std::vector<std::shared_ptr<muq::Modeling::LinearSDE>> sdes(2);
        std::vector<std::shared_ptr<muq::Utilities::LinearOperator>> linOps(2);
        Eigen::MatrixXd pinf1, pinf2;

        // Get the statespace information from each component
        std::tie(sdes.at(0), linOps.at(0), pinf1) = kernel1.GetStateSpace(sdeOptions);
        std::tie(sdes.at(1), linOps.at(1), pinf2) = kernel2.GetStateSpace(sdeOptions);

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
    LeftType  kernel1;
    RightType kernel2;
};

}
}


#endif
