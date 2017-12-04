#ifndef LINEARTRANSFORMKERNEL_H
#define LINEARTRANSFORMKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"



namespace muq
{
namespace Approximation
{
    
/**

@class LinearTransformKernel
@ingroup CovarianceKernels

Given another kernel $k_2(x,y)$ and a linear transformation $A$, this class implements a kernel of the form
\f[
k(x,y) = A * k_2(x,y) * A^T.
\f]
 */
template<typename KernelType>
class LinearTransformKernel : public KernelImpl<LinearTransformKernel<KernelType>>
{

public:

    
    LinearTransformKernel(Eigen::MatrixXd const& Ain,
		          KernelType  const& Kin) : KernelImpl<LinearTransformKernel<KernelType>>(Kin.inputDim, Ain.rows(), Kin.numParams), A(Ain), K(Kin)
    {
	assert(Ain.cols() == Kin.coDim);
    };

    virtual ~LinearTransformKernel(){};
    
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov) const
    {
	Eigen::MatrixXd tempCov(K.coDim, K.coDim);
	K.EvaluateImpl(x1,x2,tempCov);
	
	cov = A * tempCov * A.transpose();
    }

    template<typename VecType, typename MatrixType>
    inline void DerivCovarianceImpl(VecType const& x1, VecType const& x2, std::vector<unsigned> wrts, MatrixType & derivCov ) const
    {
        assert(false);
    }
        
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType &derivs) const
    {
	Eigen::MatrixXd tempDerivs(K.coDim, K.coDim);
	K.GetDerivative(x1,x2, wrt, tempDerivs);
	derivs = A * tempDerivs * A.transpose();
    }
	

    virtual Eigen::VectorXd GetParams() const override
    {
	return K.GetParams();
    }
    
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        K.SetParams(params);
    }

private:
    Eigen::MatrixXd A;
    KernelType      K;

};


template<typename KernelType>
LinearTransformKernel<KernelType> TransformKernel(Eigen::MatrixXd const& A, KernelType const& K)
{
    return LinearTransformKernel<KernelType>(A,K);
}

template<typename KernelType, typename = typename std::enable_if<std::is_base_of<KernelImpl<KernelType>, KernelType>::value, KernelType>::type>
LinearTransformKernel<KernelType> operator*(Eigen::MatrixXd const& A, KernelType const&K)
{
    return TransformKernel(A,K);
}


}
}

#endif
