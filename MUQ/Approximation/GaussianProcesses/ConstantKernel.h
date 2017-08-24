#ifndef CONSTANTKERNEL_H
#define CONSTANTKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"

#include "MUQ/Utilities/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Utilities/LinearAlgebra/ZeroOperator.h"

#include "MUQ/Modeling/LinearSDE.h"

#include <boost/property_tree/ptree.hpp>

namespace muq
{
namespace Approximation
{

/**

@class ConstantKernel
@ingroup CovarianceKernels

This class implements a constant kernel of the form
\f[
k(x,y) = \sigma^2 
\f]
where, \f$\sigma^2\f$ is the variance.

 */
class ConstantKernel : public KernelImpl<ConstantKernel>
{

public:

    ConstantKernel(unsigned              dim,
	           const double          sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : ConstantKernel(dim, sigma2In*Eigen::MatrixXd::Ones(1,1), sigmaBounds){};

    ConstantKernel(unsigned              dim,
		   std::vector<unsigned> dimInds,
	           const double          sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : ConstantKernel(dim, dimInds, sigma2In*Eigen::MatrixXd::Ones(1,1), sigmaBounds){};

    
    ConstantKernel(unsigned               dim,
	           Eigen::MatrixXd const& sigma2In,
                   const Eigen::Vector2d  sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelImpl<ConstantKernel>(dim, sigma2In.rows(), GetNumParams(sigma2In)), sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
    };

    ConstantKernel(unsigned               dim,
		   std::vector<unsigned>  dimInds,
	           Eigen::MatrixXd const& sigma2In,
                   const Eigen::Vector2d  sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelImpl<ConstantKernel>(dim, dimInds, sigma2In.rows(), GetNumParams(sigma2In)), sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
    };

    virtual ~ConstantKernel(){};
    
    template<typename VecType, typename MatrixType>
    inline void EvaluateImpl(VecType const& x1, VecType const& x2, MatrixType & cov ) const
    {
	cov = sigma2;
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	int k=0;
	for(int col=0; col<sigma2.cols(); ++col)
	{
	    for(int row=col; row<sigma2.rows(); ++row)
	    {
		if(k==wrt)
		{
		    derivs(row,col) = 1.0;
		    derivs(col,row) = 1.0;
		}
		else
		{
		    derivs(row,col) = 0.0;
		    derivs(col,row) = 0.0;
		}
	    }
	}
    }
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd output(numParams);
	int k=0;
	for(int col=0; col<sigma2.cols(); ++col)
	{
	    for(int row=col; row<sigma2.rows(); ++row)
	    {
		output(k) = sigma2(row,col);
		++k;
	    }
	}
	return output;
    }

    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        int k=0;
	for(int col=0; col<sigma2.cols(); ++col)
	{
	    for(int row=col; row<sigma2.rows(); ++row)
	    {
		sigma2(row,col) = params(k);
		++k;
	    }
	}
    }

    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions = boost::property_tree::ptree()) const override
      {
        std::shared_ptr<muq::Utilities::LinearOperator> F = std::make_shared<muq::Utilities::ZeroOperator>(sigma2.rows(),1);
        std::shared_ptr<muq::Utilities::LinearOperator> L = std::make_shared<muq::Utilities::ZeroOperator>(sigma2.rows(),1);
        Eigen::MatrixXd Q = Eigen::MatrixXd::Ones(1,1);
        
        boost::property_tree::ptree options;
        auto sde = std::make_shared<muq::Modeling::LinearSDE>(F,L,Q,options);
        
        std::shared_ptr<muq::Utilities::LinearOperator> obsOp = std::make_shared<muq::Utilities::IdentityOperator>(sigma2.rows());

        return std::make_tuple(sde, obsOp, sigma2);
      }

private:
    Eigen::MatrixXd sigma2;

    static unsigned GetNumParams(Eigen::MatrixXd const& cov)
    {
	return 0.5*cov.rows()*(cov.rows()+1);
    }
};

}
}


#endif
