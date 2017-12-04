#ifndef KERNELBASE_H
#define KERNELBASE_H


#include <assert.h>
#include <memory>
#include <iostream>
#include <fstream>

#include "MUQ/Approximation/TemplatedArrayUtilities.h"

//#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

#include "MUQ/Utilities/Exceptions.h"

#include <boost/property_tree/ptree.hpp>


namespace muq
{
namespace Modeling{
    class LinearSDE;
}
namespace Utilities{
    class LinearOperator;
}

namespace Approximation
{
    
template<typename LeftType, typename RightType>
class ProductKernel;

template<typename LeftType, typename RightType>
class SumKernel;


/** @class KernelBase
    @ingroup CovarianceKernels
    @brief Base class for all covariance kernels.
*/
class KernelBase : public std::enable_shared_from_this<muq::Approximation::KernelBase>
{

public:

    KernelBase(unsigned inputDimIn,
	       unsigned coDimIn,
	       unsigned numParamsIn) : KernelBase(inputDimIn, BuildDimInds(inputDimIn), coDimIn, numParamsIn)
    {};
    
    KernelBase(unsigned              inputDimIn,
	       std::vector<unsigned> dimIndsIn,
	       unsigned              coDimIn,
	       unsigned              numParamsIn) : dimInds(dimIndsIn), inputDim(inputDimIn), coDim(coDimIn), numParams(numParamsIn)
    {
	assert(inputDim>0);
	assert(coDim>0);
    };


    virtual ~KernelBase(){};


    //virtual std::shared_ptr<muq::Approximation::KernelBase> GetPtr() {
    //    return shared_from_this();
    //}

    /// Overridden by ProductKernel
    virtual std::vector<std::shared_ptr<KernelBase>> GetSeperableComponents() {return std::vector<std::shared_ptr<KernelBase>>(1,Clone()); };

    
    virtual Eigen::MatrixXd Evaluate(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2) const = 0;
    
    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x) const = 0;

    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x1,
                                            Eigen::MatrixXd const& x2) const = 0;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::MatrixXd             const& ys,
				Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;

    /** @brief Fills in a matrix with the first or higher order derivatives of the covariance kernel.
        @details Let \f$k(x^{(i)}),x^{(j)})\f$ be the covariance kernel taking input vectors \f$x^{(i)}\f$ and \f$x^{(j)}\f$.  This function fills in a matrix \f$\Sigma\f$ so that 
        \f[
        \Sigma_{ij} = \frac{\partial^N k(x^{(i)},x^{(j)})}{\partial x_{k_1} \partial x_{k_2} \ldots \partial x_{k_N}}, 
        \f]
where $\partial \f$x_{k_m}\f$ is the \f$k_m\f$ component of \f$x\f$.
@param[in] x1 The first location \f$x^{(i)}\f$
@param[in] x2 The second location \f$x^{(j)}\f$
@param[in] wrts  A vector containing \f$\left[k_1,k_2,\ldots, k_N\right]\f$.
@param[out] derivCov A matrix with the derivatives of the covariance.  Notice that derivCov should be sized correctly before calling this function.
    */
    virtual void FillDerivCovariance(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                     Eigen::Ref<const Eigen::VectorXd> const& x2,
                                     std::vector<unsigned>             const& wrts,
                                     Eigen::Ref<Eigen::MatrixXd>              derivCov) const = 0;
    
    virtual void FillDerivativeMatrix(Eigen::MatrixXd             const& xs,
				      unsigned                           wrt,
				      Eigen::Ref<Eigen::MatrixXd>        derivs) const = 0;
    
    virtual Eigen::MatrixXd GetDerivative(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2, int wrt) const = 0;
    
    virtual Eigen::MatrixXd GetParamBounds() const
    {
	return paramBounds;
    };


    virtual Eigen::VectorXd GetParams() const{return Eigen::VectorXd();};

    virtual void SetParams(Eigen::VectorXd const& params){};

    virtual std::shared_ptr<KernelBase> Clone() const = 0;


    /** @brief Returns a state space representation of the covariance kernel
        @details If this is a one dimensional kernel (i.e., inputDim=1 and coDim=1), this function returns a state space representation of the covariance kernel.  In particular, it returns a linear time invariant stochastic differential equation, whose solution, when started with the returned stationary covariance, provides the same information as this Gaussian process.   The first component of the vector-valued stochastic differential equation is related to the Gaussian process.  See "Kalman filtering and smoothing solutions to temporal Gaussian process regression models," by Jouni Hartikainen and Simo Sarkka, for more information.
   
    */
    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const{
        throw muq::NotImplementedError("ERROR.  The GetStateSpace() function has not been implemented in this chiled of muq::Approximation::KernelBase.");
    };
    
    
    const std::vector<unsigned> dimInds;

    const unsigned inputDim;
    const unsigned coDim;
    const unsigned numParams;
    
protected:

    Eigen::MatrixXd paramBounds;


private:
    static std::vector<unsigned> BuildDimInds(unsigned dim)
    {
	std::vector<unsigned> output(dim);
	for(int i=0; i<dim; ++i)
	    output[i] = i;

	return output;
    }

};

}
}



#endif
