#ifndef SQUAREDEXPKERNEL_H
#define SQUAREDEXPKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

    
/**

@class SquaredExpKernel
@ingroup CovarianceKernels
This class implements a kernel of the form
\f[
k(x,y) = \sigma^2 \exp\left(-\frac{1}{2}\frac{|x-y|^2}{L^2}\right)
\f]
for some variance \f$\sigma^2\f$ and lengthscale \f$L\f$.

 */
class SquaredExpKernel : public KernelImpl<SquaredExpKernel>
{

public:

    SquaredExpKernel(unsigned              dimIn,
		     std::vector<unsigned> dimInds,
		     double                sigma2In,
		     double                lengthIn,
	             Eigen::Vector2d       sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	             Eigen::Vector2d       lengthBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelImpl<SquaredExpKernel>(dimIn, dimInds, 1, 2), sigma2(sigma2In), length(lengthIn)
    {
	paramBounds.resize(2,2);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
	
	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);
    };
    
    SquaredExpKernel(unsigned        dimIn,
		     double          sigma2In,
		     double          lengthIn,
	             Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	             Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelImpl<SquaredExpKernel>(dimIn, 1, 2), sigma2(sigma2In), length(lengthIn)
    {
	paramBounds.resize(2,2);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
	
	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);
    };

    virtual ~SquaredExpKernel(){};
    
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov ) const
    {
	assert(GetShape(x2,0)==GetShape(x1,0));

	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	
	cov(0,0) = sigma2*exp(-0.5*pow(dist/length,2.0));
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	assert(wrt<numParams);

	double dist = CalcDistance( GetSlice(x1, dimInds), GetSlice(x2, dimInds) );
	
	if(wrt==0) // derivative wrt sigma2
	{
	    derivs(0,0) = exp(-0.5*pow(dist/length,2.0));
	}
	else if(wrt==1) // derivative wrt length
	{
	    derivs(0,0) = sigma2 * exp(-0.5*pow(dist/length,2.0)) * pow(dist,2.0) * pow(length, -3.0);
	}
	else
	{
	    assert(false);
	}
    }
	
    virtual Eigen::VectorXd GetParams() const override
    {
	return Eigen::Vector2d{sigma2,length};
    }
    
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        sigma2 = params(0);
	length = params(1);
    }

private:
    double sigma2;
    double length;
};

}
}


#endif
