#ifndef WHITENOISEKERNEL_H
#define WHITENOISEKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

/**

@class WhiteNoiseKernel
@ingroup CovarianceKernels

This class implements a kernel of the form
\f[
k(x,y) = \sigma^2\delta(x,y)
\f]
where \f$\delta(x,y)=1\f$ if \f$x=y\f$ and \f$0\f$ otherwise.

 */
class WhiteNoiseKernel : public KernelImpl<WhiteNoiseKernel>
{

public:

    WhiteNoiseKernel(unsigned     dim,
		     const double sigma2In,
		     const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelImpl<WhiteNoiseKernel>(dim, 1, 1), sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0) = sigmaBounds[0]; // lower bound on sigma2
	paramBounds(1) = sigmaBounds[1]; // upper bound on sigma2
    };

    virtual ~WhiteNoiseKernel(){};
    
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov) const
    {
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	cov(0,0) = (dist<1e-14) ? sigma2 : 0.0;
    }

    template<typename VecType, typename MatrixType>
    inline void DerivCovarianceImpl(VecType const& x1, VecType const& x2, std::vector<unsigned> wrts, MatrixType & derivCov ) const
    {
        assert(false);
    }
        
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType &derivs) const
    {
	assert(wrt==0);
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));

	derivs(0,0) = (dist<1e-14) ? 1.0 : 0.0;
    }
	

    virtual Eigen::VectorXd GetParams() const override
    {
	return sigma2*Eigen::VectorXd::Ones(1);
    }
    
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        sigma2 = params(0);
    }

private:
    double sigma2;

};


}
}

#endif
