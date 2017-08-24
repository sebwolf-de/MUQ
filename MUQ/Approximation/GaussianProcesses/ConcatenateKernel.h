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
    class ConcatenateKernel : public KernelBase
    {

    public:
        ConcatenateKernel(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn);
        virtual ~ConcatenateKernel(){};
    
        virtual Eigen::MatrixXd Evaluate(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2) const override;
        
        virtual void FillDerivativeMatrix(Eigen::MatrixXd             const& xs,
                                          unsigned                           wrt,
                                          Eigen::Ref<Eigen::MatrixXd>        derivs) const override;
        
        
        virtual void FillCovariance(Eigen::MatrixXd             const& xs,
                                    Eigen::MatrixXd             const& ys,
                                    Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;
        
        virtual void FillCovariance(Eigen::MatrixXd             const& xs,
                                    Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;
        
    
        virtual Eigen::MatrixXd GetParamBounds() const override
        {
            Eigen::MatrixXd bounds(2,this->numParams);

            int currCol = 0;
            for(int i=0; i<kernels.size(); ++i){
              bounds.block(0,currCol, 2, kernels.at(i)->numParams) = kernels.at(i)->GetParamBounds();
              currCol += kernels.at(i)->numParams;
            }

            return bounds;
        };
	
        virtual Eigen::VectorXd GetParams() const override
        {
            Eigen::VectorXd params(this->numParams);

            int currCol = 0;
            for(int i=0; i<kernels.size(); ++i){
              params.segment(currCol,kernels.at(i)->numParams) = kernels.at(i)->GetParams();
              currCol += kernels.at(i)->numParams;
            }

            return params;
        };
        
        virtual void SetParams(Eigen::VectorXd const& params) override
        {
            int currCol = 0;
            for(int i=0; i<kernels.size(); ++i){
              kernels.at(i)->SetParams(params.segment(currCol,kernels.at(i)->numParams));
              currCol += kernels.at(i)->numParams;
            }
        }

    
    private:
        static int GetNumCodim(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn);
        static int GetNumParams(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn);

        std::vector<std::shared_ptr<KernelBase>> kernels;
    };


    
template<class KernelType1, class KernelType2>
ConcatenateKernel Concatenate(KernelType1 const& kernel1, KernelType2 const& kernel2)
{
    std::vector<std::shared_ptr<KernelBase>> kernels(2);
    kernels.at(0) = kernel1.Clone();
    kernels.at(1) = kernel2.Clone();

    return ConcatenateKernel(kernels);
}




} // namespace Approximation
} // namespace muq


#endif
