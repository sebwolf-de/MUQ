#ifndef PERIODICKERNEL_H
#define PERIODICKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

    /** 

@class PeriodicKernel

\f[
k(x,y) = \exp\left[-\frac{2}{L^2} \sin^2\left(\frac{\pi}{p}(x-y)\right)\right]

 */
class PeriodicKernel : public KernelImpl<PeriodicKernel>
{

public:

     PeriodicKernel(unsigned dim,
		    std::vector<unsigned> dimInds,
		    const double sigma2In,
		    const double lengthIn,
		    const double periodIn,
                    const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	            const Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()},
		    const Eigen::Vector2d periodBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelImpl<PeriodicKernel>(dim, dimInds, 1 , 3),
	                                                                                                     sigma2(sigma2In),
	                                                                                                     length(lengthIn),
	                                                                                                     period(periodIn)
    {
	SetupBounds(sigmaBounds, lengthBounds, periodBounds);
    };

    virtual ~PeriodicKernel(){};
    
    PeriodicKernel(unsigned dim,
		   const double sigma2In,
		   const double lengthIn,
		   const double periodIn,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	           const Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()},
	           const Eigen::Vector2d periodBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelImpl<PeriodicKernel>(dim, 1 , 3),
	                                                                                                    sigma2(sigma2In),
	                                                                                                    length(lengthIn),
	                                                                                                    period(periodIn)
    {
	SetupBounds(sigmaBounds, lengthBounds, periodBounds);
    };

    
    template<typename VecType, typename MatType>
    inline void EvaluateImpl(VecType const& x1, VecType const& x2, MatType & cov) const
    {
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	cov(0,0) = sigma2 * exp(-2.0 * pow(sin(pi*dist/period),2.0) / pow(length,2.0));
    }


    template<typename VecType1, typename VecType2, typename MatType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < this->numParams);

	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	
	if( wrt == 0)
	{
	    derivs(0,0) = exp(-2.0 * pow(sin(pi*dist/period),2.0) / pow(length,2.0));

	}
	else if(wrt==1)
	{
	    derivs(0,0) = 4.0 * sigma2 * exp(-2.0 * pow(sin(pi*dist/period),2.0) / pow(length,2.0)) *  std::pow(sin(pi*dist/period),2.0) * std::pow(length,-3.0);
	}
	else if(wrt==2)
	{
	    double temp = pi*dist/period;
	    derivs(0,0) = 4.0 * sigma2 * exp(-2.0 * pow(sin(temp),2.0) / pow(length,2.0)) * sin(temp) * cos(temp) * pi * dist * std::pow(length*period,-2.0);
	}
	else
	{
	    assert(false);
	}
    }
	
    virtual Eigen::VectorXd GetParams() const override
    {
	return Eigen::Vector3d{sigma2, length, period};
    };
    
    virtual void  SetParams(Eigen::VectorXd const& params) override
    {
        sigma2 = params(0);
	length = params(1);
	period = params(2);
    }

    
private:
    double sigma2;
    double length;
    double period;

    const double pi = 4.0 * atan(1.0); //boost::math::constants::pi<double>();


    void SetupBounds(Eigen::Vector2d const& sigmaBounds,
		     Eigen::Vector2d const& lengthBounds,
		     Eigen::Vector2d const& periodBounds)
    {
	paramBounds.resize(2,3);

	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);

	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);

	paramBounds(0,2) = periodBounds(0);
	paramBounds(1,2) = periodBounds(1);
    };
};
    
}
}


#endif
