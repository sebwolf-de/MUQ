#ifndef CONCATENATEKERNEL_H
#define CONCATENATEKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


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
    class ConcatenateKernel : public KernelBase
    {

    public:
        ConcatenateKernel(std::shared_ptr<KernelBase> const& kernel1In,
                          std::shared_ptr<KernelBase> const& kernel2In) : KernelBase(kernel1In->inputDim,
                                                                                     kernel1In->coDim + kernel2In->coDim,
                                                                                     kernel1In->numParams + kernel2In->numParams),
                                                                          kernel1(kernel1In),
                                                                          kernel2(kernel2In)
        {
            assert(kernel1->inputDim == kernel2->inputDim);
        };

        virtual ~ConcatenateKernel(){};

        virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                               Eigen::Ref<const Eigen::VectorXd> const& x2,
                               Eigen::Ref<const Eigen::VectorXd> const& params,
                               Eigen::Ref<Eigen::MatrixXd>              block) const override
        {
          block = Eigen::MatrixXd::Zero(coDim, coDim);
          kernel1->FillBlock(x1, x2, params, block.block(0,0,kernel1->coDim, kernel1->coDim));
          kernel2->FillBlock(x1, x2, params, block.block(kernel1->coDim,kernel1->coDim,kernel2->coDim, kernel2->coDim));
        }


        // template<typename VecType1, typename VecType2, typename MatType>
        //     inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
        // {
        //     assert(wrt < this->numParams);
        //
        //     // Initialize the derivative matrix to 0
        //     for(int j=0; j<derivs.cols(); ++j)
        //     {
        //         for(int i=0; i<derivs.rows(); ++i)
        //             derivs(i,j) = 0.0;
        //     }
        //
        //     if(wrt < kernel1.numParams )
        //     {
        //         auto block = GetBlock(derivs, 0, 0, kernel1.coDim, kernel1.coDim);
        //         return kernel1.GetDerivative(x1, x2, wrt, block);
        //     }
        //     else
        //     {
        //         auto block = GetBlock(derivs, kernel1.coDim, kernel1.coDim, kernel2.coDim, kernel2.coDim);
        //         return kernel2.GetDerivative(x1, x2, wrt-kernel1.numParams, block);
        //     }
        // }

    private:
        std::shared_ptr<KernelBase> kernel1;
        std::shared_ptr<KernelBase> kernel2;
    };



template<class KernelType1, class KernelType2>
ConcatenateKernel<KernelType1, KernelType2> Concatenate(KernelType1 const& kernel1, KernelType2 const& kernel2)
{
    return ConcatenateKernel(kernel1.Clone(), kernel2.Clone());
}




} // namespace Approximation
} // namespace muq


#endif
