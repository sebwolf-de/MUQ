#ifndef CONCATENATEKERNEL_H
#define CONCATENATEKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

    /**
       @class ConcatenateKernel
       @ingroup CovarianceKernels
       @brief Concatenates two kernels together.
       @details Let \f$k_1(x,x^\prime)\f$ and \f$k_2(x,x^\prime)\f$ be two different covariance kernels with the same inputs.  This class describes a concatenated kernel of the form 
\f[
k(x,x^\prime) = \left[\begin{array}{cc}k_1(x,x^\prime) & 0\\ 0 & k_2(x,x^\prime)\end{array}\right].
\f]
    */
    template<typename KernelType1, typename KernelType2>
    class ConcatenateKernel : public KernelImpl<ConcatenateKernel<KernelType1,KernelType2>>
    {

    public:
        ConcatenateKernel(KernelType1 const& kernel1In, KernelType2 const& kernel2In) : KernelImpl<ConcatenateKernel<KernelType1,KernelType2>>(kernel1In.inputDim,
                                                                                                                                               kernel1In.coDim + kernel2In.coDim,
                                                                                                                                               kernel1In.numParams + kernel2In.numParams),
            kernel1(kernel1In),
            kernel2(kernel2In)
        {
            assert(kernel1.inputDim == kernel2.inputDim);
        };

        virtual ~ConcatenateKernel(){};
    
        template<typename VecType, typename MatType>
            inline void EvaluateImpl(VecType const& x1, VecType  const& x2 , MatType & cov) const
        {
            auto block1 = GetBlock(cov, 0, 0, kernel1.coDim, kernel1.coDim);
            auto block2 = GetBlock(cov, kernel1.coDim, kernel1.coDim, kernel2.coDim, kernel2.coDim);
        
            kernel1.EvaluateImpl(x1, x2, block1);
            kernel2.EvaluateImpl(x1, x2, block2);
        }

        template<typename VecType1, typename VecType2, typename MatType>
            inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
        {
            assert(wrt < this->numParams);

            // Initialize the derivative matrix to 0
            for(int j=0; j<derivs.cols(); ++j)
            {
                for(int i=0; i<derivs.rows(); ++i)
                    derivs(i,j) = 0.0;
            }
        
            if(wrt < kernel1.numParams )
            {
                auto block = GetBlock(derivs, 0, 0, kernel1.coDim, kernel1.coDim);
                return kernel1.GetDerivative(x1, x2, wrt, block);
            }
            else
            {
                auto block = GetBlock(derivs, kernel1.coDim, kernel1.coDim, kernel2.coDim, kernel2.coDim);
                return kernel2.GetDerivative(x1, x2, wrt-kernel1.numParams, block);
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

    
    private:
        KernelType1  kernel1;
        KernelType2 kernel2;
    };


    
template<class KernelType1, class KernelType2>
ConcatenateKernel<KernelType1, KernelType2> Concatenate(KernelType1 const& kernel1, KernelType2 const& kernel2)
{
    return ConcatenateKernel<KernelType1, KernelType2>(kernel1, kernel2);
}




} // namespace Approximation
} // namespace muq


#endif
