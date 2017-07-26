#ifndef PRODUCTKERNEL_H
#define PRODUCTKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"

#include "MUQ/Utilities/LinearAlgebra/KroneckerProductOperator.h"
#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Utilities/LinearAlgebra/BlockRowOperator.h"
#include "MUQ/Utilities/LinearAlgebra/SumOperator.h"
#include "MUQ/Utilities/LinearAlgebra/IdentityOperator.h"

#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

#include "MUQ/Utilities/Exceptions.h"

namespace muq
{
namespace Approximation
{

template<typename LeftType, typename RightType>
std::shared_ptr<StateSpaceGP> GetProductStateSpace(LeftType const& kernel1, RightType const& kernel2, boost::property_tree::ptree sdeOptions)
{
    throw muq::NotImplementedError("ERROR.  The GetStateSpace() function has not been implemented in this child of muq::Approximation::KernelBase.");
}


// See "Explicit Link Between Periodic 
template<typename RightType>
std::shared_ptr<StateSpaceGP> GetProductStateSpace(PeriodicKernel const& kernel1, RightType const& kernel2, boost::property_tree::ptree sdeOptions)
{

    std::shared_ptr<StateSpaceGP> periodicGP = kernel1.GetStateSpace(sdeOptions);
    auto periodicSDE = periodicGP->GetSDE();

    auto periodicF = std::dynamic_pointer_cast<muq::Utilities::BlockDiagonalOperator>(periodicSDE->GetF());
    assert(periodicF);
    
    auto periodicL = std::dynamic_pointer_cast<muq::Utilities::BlockDiagonalOperator>(periodicSDE->GetL());
    assert(periodicL);
    
    std::shared_ptr<StateSpaceGP> otherGP = kernel2.GetStateSpace(sdeOptions);
    auto otherSDE = otherGP->GetSDE();
    auto otherF = otherSDE->GetF();
    auto otherL = otherSDE->GetL();
    auto otherH = otherGP->GetObs();
    
    /// Construct the new F operator
    std::vector<std::shared_ptr<muq::Utilities::LinearOperator>> newBlocks( periodicF->GetBlocks().size() );
    for(int i=0; i<newBlocks.size(); ++i)
        newBlocks.at(i) = muq::Utilities::KroneckerSum(otherF, periodicF->GetBlock(i) );

    auto newF = std::make_shared<muq::Utilities::BlockDiagonalOperator>(newBlocks);
    
    /// Construct the new L operator
    for(int i=0; i<newBlocks.size(); ++i)
        newBlocks.at(i) = std::make_shared<muq::Utilities::KroneckerProductOperator>(otherL, periodicL->GetBlock(i) );

    auto newL = std::make_shared<muq::Utilities::BlockDiagonalOperator>(newBlocks);

    /// Construct the new H operator
    Eigen::MatrixXd Hblock(1,2);
    Hblock << 1.0, 0.0;
    
    for(int i=0; i<newBlocks.size(); ++i)
        newBlocks.at(i) = std::make_shared<muq::Utilities::KroneckerProductOperator>(otherH, muq::Utilities::LinearOperator::Create(Hblock) );

    auto newH = std::make_shared<muq::Utilities::BlockRowOperator>(newBlocks);

    // Construct Pinf
    Eigen::MatrixXd periodicP = periodicGP->GetCov();
    Eigen::MatrixXd otherP = otherGP->GetCov();
    
    Eigen::MatrixXd Pinf = Eigen::MatrixXd::Zero(periodicP.rows()*otherP.rows(), periodicP.cols()*otherP.cols());
    for(int i=0; i<newBlocks.size(); ++i)
        Pinf.block(2*i*otherP.rows(), 2*i*otherP.rows(), 2*otherP.rows(), 2*otherP.cols()) = muq::Utilities::KroneckerProduct(otherP, periodicP.block(2*i,2*i,2,2));

    // Construct Q
    Eigen::MatrixXd const& otherQ = otherSDE->GetQ(); 
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(otherQ.rows()*periodicP.rows(), otherQ.cols()*periodicP.cols());
    for(int i=0; i<newBlocks.size(); ++i)
        Q.block(2*i*otherQ.rows(), 2*i*otherQ.rows(), 2*otherQ.rows(), 2*otherQ.cols()) = muq::Utilities::KroneckerProduct(otherQ, periodicP.block(2*i,2*i,2,2));

    // Construct the new statespace GP
    auto newSDE = std::make_shared<muq::Modeling::LinearSDE>(newF, newL, Q, sdeOptions);
    auto newGP = std::make_shared<StateSpaceGP>(newSDE, newH, Pinf);
    
    return newGP;
}

template<typename LeftType>
std::shared_ptr<StateSpaceGP> GetProductStateSpace(LeftType const& kernel1, PeriodicKernel const& kernel2, boost::property_tree::ptree sdeOptions)
{
    return GetProductStateSpace(kernel2, kernel1, sdeOptions);
}



    
/**

@class ProductKernel
#ingroup CovarianceKernels
\f[
k(x,y) = k_1(x,y)*k_2(x,y)
\f]

 */
template<typename LeftType, typename RightType>
class ProductKernel : public KernelImpl<ProductKernel<LeftType,RightType>>
{

public:
    ProductKernel(LeftType const& kernel1In, RightType const& kernel2In) : KernelImpl<ProductKernel<LeftType,RightType>>(kernel1In.inputDim,
														     std::max(kernel1In.coDim, kernel2In.coDim),
														     kernel1In.numParams + kernel2In.numParams),
	                                                               kernel1(kernel1In),
	                                                               kernel2(kernel2In)
    {
	assert((kernel1.coDim==kernel2.coDim) | (kernel1.coDim==1) | (kernel2.coDim==1));
    };

    virtual ~ProductKernel(){};
    
    template<typename VecType, typename MatType>
    inline void EvaluateImpl(VecType const& x1, VecType const& x2, MatType & cov ) const
    {
	Eigen::MatrixXd temp1(kernel1.coDim, kernel1.coDim);
	Eigen::MatrixXd temp2(kernel2.coDim, kernel2.coDim);

	kernel1.EvaluateImpl(x1,x2, temp1);
	kernel2.EvaluateImpl(x1,x2, temp2);

	if(kernel1.coDim==kernel2.coDim)
	{
	    cov = Eigen::MatrixXd(temp1.array() * temp2.array());
	}
	else if(kernel1.coDim==1)
	{
	    cov = temp1(0,0)*temp2;
	}
	else if(kernel2.coDim==1)
	{
	    cov = temp2(0,0)*temp1;
	}
	else
	{
	    std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
	    assert(false);
	}
	
    }

    template<typename VecType1, typename VecType2, typename MatType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < this->numParams);

	Eigen::MatrixXd temp1(kernel1.coDim, kernel1.coDim);
	Eigen::MatrixXd temp2(kernel2.coDim, kernel2.coDim);

	if(wrt < kernel1.numParams )
	{	    
	    kernel1.GetDerivative(x1, x2, wrt, temp1);
	    kernel2.EvaluateImpl(x1,x2, temp2);
	}
	else
	{
	    kernel1.EvaluateImpl(x1,x2, temp1);
	    kernel2.GetDerivative(x1,x2, wrt-kernel1.numParams, temp2);
	}

	if(kernel1.coDim==kernel2.coDim)
	{
	    derivs = Eigen::MatrixXd(temp1.array() * temp2.array());
	}
	else if(kernel1.coDim==1)
	{
	    derivs = temp1(0,0) * temp2;
	}
	else if(kernel2.coDim==1)
	{
            derivs = temp2(0,0) * temp1;
	}
	else
	{
	    std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
	    assert(false);
	}
	
    }

    virtual std::vector<std::shared_ptr<KernelBase>> GetSeperableComponents() override
    {
        // Check if the dimensions of the components are distinct
        bool isSeperable = true;
        for(unsigned leftDim : kernel1.dimInds)
        {
            for(unsigned rightDim : kernel2.dimInds)
            {
                if(leftDim==rightDim)
                {
                    isSeperable = false;
                    break;
                }
            }
            
            if(!isSeperable)
                break;
        }

        if(isSeperable)
        {
            std::vector<std::shared_ptr<KernelBase>> output, output2;
            output = kernel1.GetSeperableComponents();
            output2 = kernel2.GetSeperableComponents();
            output.insert(output.end(), output2.begin(), output2.end());

            return output;
        }
        else
        {
            return std::vector<std::shared_ptr<KernelBase>>(1, KernelBase::GetPtr());
        }
            
    };
    
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

    virtual std::shared_ptr<StateSpaceGP> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const override{
        return GetProductStateSpace(kernel1, kernel2, sdeOptions);
    };
    
protected:
    LeftType  kernel1;
    RightType kernel2;

};


}
}


#endif
