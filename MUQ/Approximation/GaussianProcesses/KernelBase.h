#ifndef KERNELBASE_H
#define KERNELBASE_H


#include <assert.h>
#include <memory>
#include <iostream>
#include <fstream>

#include "MUQ/Approximation/GaussianProcesses/Utilities.h"

namespace muq
{
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
class KernelBase
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
    
    virtual Eigen::MatrixXd Evaluate(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2) const = 0;
    
    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x) const = 0;

    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x1,
                                            Eigen::MatrixXd const& x2) const = 0;
    
    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::MatrixXd             const& ys,
				Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;

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
