#ifndef COREGIONALKERNEL_H
#define COREGIONALKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

/** 

@class CoregionalKernel
@ingroup CovarianceKernels

This kernel supports coregionalization for modeling vector-valued Gaussian processes.

 */
class CoregionalKernel : public KernelImpl<CoregionalKernel>
{

public:
    
    CoregionalKernel(unsigned               dim,
		     Eigen::MatrixXd const& Gamma,
		     std::vector<std::shared_ptr<KernelBase>> const& kernelsIn) : KernelImpl<CoregionalKernel>(dim, Gamma.rows(), GetNumParams(kernelsIn)),
                                                                                  kernels(kernelsIn)
    {
	// Make sure the matrix and kernels are the same size
	assert(Gamma.cols()==kernelsIn.size());

	// Compute the eigenvalue decomposition of the covariance Gamma to get the matrix A
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver;
	eigSolver.compute(Gamma, Eigen::ComputeEigenvectors);
	A = eigSolver.eigenvectors()*eigSolver.eigenvalues().cwiseSqrt().asDiagonal();
    };

    virtual ~CoregionalKernel(){};

    template<typename VecType, typename MatrixType>
    inline void EvaluateImpl(VecType const& x1, VecType const& x2, MatrixType & cov ) const
    {
	Eigen::VectorXd rhoVec(this->coDim);
        for(int i=0; i<this->coDim; ++i)
            rhoVec(i) = kernels.at(i)->Evaluate(x1,x2)(0,0);
        
	cov = ( A * rhoVec.asDiagonal() * A.transpose() ).eval();
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	unsigned cumParams = 0;
        unsigned kernelInd = 0;
        Eigen::MatrixXd kernelDerivs;
        
        for(unsigned i=0; i<this->coDim; ++i)
        {
            if(wrt < cumParams + kernels.at(i)->numParams)
            {
                kernelDerivs = kernels.at(i)->GetDerivative(x1, x2, wrt-cumParams);
                kernelInd = i;
                break;
            }
            else
            {
                cumParams += kernels.at(i)->numParams;
            }
        }

        assert(kernelDerivs.rows()==1);
        assert(kernelDerivs.cols()==1);
        
	derivs = Eigen::MatrixXd(A.col(kernelInd) * kernelDerivs(0,0) * A.col(kernelInd).transpose());
    }


    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(this->numParams);
        unsigned cumParams = 0;
        for(int i=0; i<this->coDim; ++i)
        {
            params.segment(cumParams, kernels.at(i)->numParams) = kernels.at(i)->GetParams();
            cumParams += kernels.at(i)->numParams;
        }
        
	return params;
    }

    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        unsigned cumParams = 0;
        for(int i=0; i<this->coDim; ++i)
        {
            kernels.at(i)->SetParams( params.segment(cumParams, kernels.at(i)->numParams));
            cumParams += kernels.at(i)->numParams;
        }
        
    }
	
private:

    Eigen::MatrixXd A;
    std::vector<std::shared_ptr<KernelBase>> kernels;

    static unsigned GetNumParams(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn)
    {
        unsigned allParams = 0;
        for(unsigned i = 0; i<kernelsIn.size(); ++i)
            allParams += kernelsIn.at(i)->numParams;
        return allParams;
    };

};



template<class KernelType1, class... KTypes>
void FillKernelVector(std::vector<std::shared_ptr<KernelBase>> & vec, const KernelType1& kernel1, const KTypes&... kernels)
{
    vec.push_back(kernel1.Clone());
    FillKernelVector(vec, kernels...);
}

template<class KernelType1>
void FillKernelVector(std::vector<std::shared_ptr<KernelBase>> & vec, const KernelType1& kernel1)
{
    KernelBase const& baseRef = kernel1;
    vec.push_back(baseRef.Clone());
}


/**

Create a coregional kernel
@param[in] coCov The joint covariance of the vector-valued process
@param[in] kernels Covariance kernels for the principal components of the covariance.
 
*/
template<class KernelType1, class KernelType2, typename = typename std::enable_if<std::is_base_of<KernelBase, KernelType2>::value>, class... KTypes>
CoregionalKernel CoregionTie(Eigen::MatrixXd const& coCov, const KernelType1& kernel1, const KernelType2& kernel2, const KTypes&... kernels)
{
    std::vector<std::shared_ptr<KernelBase>> vec;
    FillKernelVector(vec, kernel1, kernel2, kernels...);
    
    return CoregionalKernel(kernel1.inputDim, coCov, vec);

}
template<class KernelType1>
CoregionalKernel CoregionTie(Eigen::MatrixXd const& coCov, const KernelType1& kernel1, int numRepeat)
{
    std::vector<std::shared_ptr<KernelBase>> vec(numRepeat);
    KernelBase const& baseRef = kernel1;
    for(int i=0; i<numRepeat; ++i)
        vec.at(i) = baseRef.Clone();

    return CoregionalKernel(kernel1.inputDim, coCov, vec);
}



}
}


#endif
