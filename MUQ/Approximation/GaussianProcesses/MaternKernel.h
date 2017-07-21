#ifndef MATERNKERNEL_H
#define MATERNKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"

#include <cmath>
#include <stdexcept>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace muq
{
namespace Approximation
{

    
/**

@class MaternKernel
@ingroup CovarianceKernels
This class implements a kernel of the form
\f[
k(x,y) = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)}\left(\frac{\sqrt{2\nu}}{l}\tau\right)^\nu K_\nu\left(\frac{\sqrt{2\nu}}{l}\tau\right),
\f]
where \f$\Gamma(\cdot)\f$ is the gamma function, \f$K_\nu(\cdot)\f$ is a modified Bessel function of the second kind, \f$\nu\f$ is a smoothness parameter, \f$l\f$ is a lengthscale, and \f$\sigma^2\f$ is the variance.  Note that we only allow values of \f$\nu\f$ that take the form \f$i-0.5\f$ for some positive integer \f$i\f$.  Typical choices are \f$\nu=1/2\f$, \f$\nu=3/2\f$, and \f$\nu=5/2\f$.   In the limit, \f$\nu\rightarrow\infty\f$, this Matern kernel converges to the squared exponential kernel (muq::Approximation::SquaredExpKernel) and for \f$\nu=1/2\f$, the Matern kernel is equivalent to the exponential kernel associated with the Ornstein-Uhlenbeck process. 

Note: the smoothness parameter \f$\nu\f$ is not optimized as a hyperparameter.
 */
class MaternKernel : public KernelImpl<MaternKernel>
{

public:

    MaternKernel(unsigned              dimIn,
                 std::vector<unsigned> dimInds,
                 double                sigma2In,
                 double                lengthIn,
                 double                nuIn,
                 Eigen::Vector2d       sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
                 Eigen::Vector2d       lengthBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelImpl<MaternKernel>(dimIn, dimInds, 1, 2), sigma2(sigma2In), length(lengthIn), nu(nuIn), scale(std::pow(2.0, 1.0-nuIn)/boost::math::tgamma(nuIn))
    {
        CheckNu();
        
	paramBounds.resize(2,2);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
	
	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);
    };
    
    MaternKernel(unsigned        dimIn,
                 double          sigma2In,
                 double          lengthIn,
                 double          nuIn,
                 Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
                 Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelImpl<MaternKernel>(dimIn, 1, 2), sigma2(sigma2In), length(lengthIn), nu(nuIn), scale(std::pow(2.0, 1.0-nuIn)/boost::math::tgamma(nuIn))
    {
        CheckNu();
        
	paramBounds.resize(2,2);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
	
	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);
    };

    virtual ~MaternKernel(){};
    
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov ) const
    {
	assert(GetShape(x2,0)==GetShape(x1,0));

	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));

        double temp = sqrt(2.0*nu)*dist/length;
	cov(0,0) = sigma2 * scale * std::pow(temp, nu) * boost::math::cyl_bessel_k(nu, temp);
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	assert(wrt<numParams);

	double dist = CalcDistance( GetSlice(x1, dimInds), GetSlice(x2, dimInds) );
	double temp = sqrt(2.0*nu)*dist/length;
        
	if(wrt==0) // derivative wrt sigma2
	{
	    derivs(0,0) = scale * std::pow(temp, nu) * boost::math::cyl_bessel_k(nu, temp);
	}
	else if(wrt==1) // derivative wrt length
	{
            double dtemp_dlength = -1.0*sqrt(2.0*nu)*dist/(length*length);
            double part1 = nu*std::pow(temp,nu-1)*dtemp_dlength;
            double part2 = -0.5*(boost::math::cyl_bessel_k(nu-1.0, temp) + boost::math::cyl_bessel_k(nu+1.0,temp) )*dtemp_dlength;
	    derivs(0,0) = sigma2 * scale * (part1*boost::math::cyl_bessel_k(nu, temp) + part2*std::pow(temp, nu));
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
    double nu;
    double scale;

    void CheckNu() const{

        if(nu<=0)
            throw std::invalid_argument("The value of nu must be greater than 0.");

        if((std::round(nu-0.5)-(nu-0.5)) > 4.0*std::numeric_limits<double>::epsilon())
            throw std::invalid_argument("The value of nu must take the form nu=i-0.5 for a positive integer i.");

    };
};

}
}


#endif
