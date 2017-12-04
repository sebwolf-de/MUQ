#ifndef MATERNKERNEL_H
#define MATERNKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"

#include <cmath>
#include <stdexcept>

#include <boost/property_tree/ptree_fwd.hpp>


#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/constants/constants.hpp>



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
                 Eigen::Vector2d       lengthBounds = {1e-10, std::numeric_limits<double>::infinity()});
    
    
    MaternKernel(unsigned        dimIn,
                 double          sigma2In,
                 double          lengthIn,
                 double          nuIn,
                 Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
                 Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()});

    
    virtual ~MaternKernel(){};

    
    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov ) const
    {
	assert(GetShape(x2,0)==GetShape(x1,0));

	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));

        if(dist < 4.0*std::numeric_limits<double>::epsilon()){
            cov(0,0) = sigma2;
        }else{
        
            double temp = sqrt(2.0*nu)*dist/length;
            cov(0,0) = sigma2 * scale * std::pow(temp, nu) * boost::math::cyl_bessel_k(nu, temp);
        }
    }

    
    template<typename VecType, typename MatrixType>
    inline void DerivCovarianceImpl(VecType const& x1, VecType const& x2, std::vector<unsigned> wrts, MatrixType & derivCov ) const
    {
        assert(false);
        
        //assert(GetShape(x2,0)==GetShape(x1,0));
        
        // Check to make sure all of the wrts are in dimInds.  If not, we the derivative will be zero.
        
        //double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));

        // get the derivative of the distance wrt to the variable

        /* // get the derivative of the kernel wrt to distance */
        /* const unsigned derivOrder = wrts.size(); */

        /* double x = sqrt(2*nu)*dist/length; */

        /* double bessel0 = boost::math::cyl_bessel_k(nu, x); */
        /* double bessel1 = boost::math::cyl_bessel_k(nu-1, x); */
        
        /* double c = std::pow(x, nu) * bessel1; */

        /* //double dc_dx = - std::pow<double>(x,nu) * bessel1; */
        /* double d2c_dx2 = (1.0-2.0*nu)*std::pow<double>(x,nu-1)*bessel1 + std::pow<double>(x,nu)*bessel0; */
        /* double d3c_dx3 = (2.0*nu-1.0) * std::pow<double>(x, nu-1)*bessel0 - std::pow<double>(x,nu-2.0)*(4.0*nu*nu-6.0*nu+x*x+2.0)*bessel1; */
        /* double d4c_dx4 = std::pow<double>(x,nu-2.0)*(4.0*nu*nu-8.0*nu+x*x+3)*bessel0 - 2.0*(2.0*nu-1.0)*std::pow<double>(x,nu-3.0)*(2.0*nu*nu-5.0*nu+x*x+3)*bessel1; */

        /* double dx_dd = sqrt(2.0*nu)/length; */
        // note that higher order derivatives of x wrt d are zero

        /* double dc_dd = dc_dx * dx_dd; */
        /* double d2c_dd2 = d2c_dx2 * std::pow<double>(dx_dd,2.0); */
        /* double d3c_dd3 = d3c_dx3 * std::pow<double>(dx_dd,3.0); */
        /* double d3c_dd4 = d4c_dx4 * std::pow<double>(dx_dd,4.0); */

        // Now that we have the derivatives wrt the distance, compute the derivatives wrt the inputs
        
        
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

    
    virtual Eigen::VectorXd GetParams() const override;
    
    virtual void SetParams(Eigen::VectorXd const& params) override;

    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions = boost::property_tree::ptree()) const override;

private:
    
    double sigma2;
    double length;
    double nu;
    double scale;

    void CheckNu() const;
};

}
}


#endif
