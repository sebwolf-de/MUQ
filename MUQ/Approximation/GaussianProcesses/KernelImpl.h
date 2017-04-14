#ifndef KERNELIMPL_H
#define KERNELIMPL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"



namespace muq
{
namespace Approximation
{

/**

\class KernelImpl
\ingroup CovarianceKernels
\brief Base class in CRTP pattern for covariance kernels
\details This class provides common functionality (such as computing Covariance matrices) for all covariance kernels.  It uses the curiously recurring template pattern and requires that child classes implement the following functions
- void EvaluateImpl(VectorType1, VectorType2, MatType)
- void GetDerivative(VectorType1, VectorType2, MatType)
- Eigen::VectorXd GetParams()
- void SetParams(Eigen::VectorXd)
*/
template<typename ChildType>
class KernelImpl : public KernelBase
{

    
public:

    KernelImpl(unsigned inputDimIn,
	       unsigned coDimIn,
	       unsigned numParamsIn) :KernelBase(inputDimIn, coDimIn, numParamsIn){};
    
    KernelImpl(unsigned              inputDimIn,
	       std::vector<unsigned> dimIndsIn,
	       unsigned              coDimIn,
	       unsigned              numParamsIn) : KernelBase(inputDimIn, dimIndsIn, coDimIn, numParamsIn){};


    virtual ~KernelImpl(){};
    
    virtual std::shared_ptr<KernelBase> Clone() const override
    {
	return std::make_shared<ChildType>(static_cast<ChildType const &>(*this));
    }

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::MatrixXd             const& ys,
				Eigen::Ref<Eigen::MatrixXd>        cov) const override
    {
	FillCovarianceImpl(xs,ys,cov);
    };

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::Ref<Eigen::MatrixXd>        cov) const override
    {
	FillCovarianceImpl(xs,cov);
    };

    virtual void FillDerivativeMatrix(Eigen::MatrixXd             const& xs,
				      unsigned                           wrt,
				      Eigen::Ref<Eigen::MatrixXd>        derivs) const override
    {
	FillDerivativeMatrixImpl(xs,wrt, derivs);
    };

    
    virtual Eigen::MatrixXd GetDerivative(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2, int wrt) const override
    {
        Eigen::MatrixXd output(coDim, coDim);
        reinterpret_cast<const ChildType*>(this)->GetDerivative(x1,x2,wrt,output);
        return output;
    };
    
    template<typename RightType>
    ProductKernel<ChildType, RightType> operator*(RightType const& kernel2) const
    {
	return ProductKernel<ChildType,RightType>(*reinterpret_cast<const ChildType*>(this),kernel2);
    }

    template<typename RightType>
    SumKernel<ChildType,RightType> operator+(RightType const& kernel2) const
    {
	return SumKernel<ChildType,RightType>(*reinterpret_cast<const ChildType*>(this),kernel2);
    }

    virtual Eigen::MatrixXd Evaluate(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2) const override
    {

        Eigen::MatrixXd output(coDim, coDim);
        static_cast<const ChildType*>(this)->EvaluateImpl( x1, x2, output);
        return output;
    };
    
    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x) const override
    {
	Eigen::MatrixXd output(coDim*x.cols(), coDim*x.cols());
        FillCovariance(x, output);
        return output;
    }


    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x1,
                                            Eigen::MatrixXd const& x2) const override
    {
	Eigen::MatrixXd output(coDim*x1.cols(), coDim*x2.cols());
        FillCovariance(x1, x2, output);
	return output;
    }

    
    template<typename PosMatrixType, typename DerivMatrixType>
    void FillDerivativeMatrixImpl(PosMatrixType const& xs, int wrt, DerivMatrixType &derivs) const
    {
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

    	assert(GetShape(derivs, 0) == coDim*GetShape(xs,1));
    	assert(GetShape(derivs, 1) == coDim*GetShape(xs,1));

	for(unsigned col=0; col<N1; ++col)
	{
	    for(unsigned row=col; row<N1; ++row)
	    {
		auto block = GetBlock(derivs, row*coDim, col*coDim, coDim, coDim);
    		static_cast<const ChildType*>(this)->GetDerivative( GetColumn(xs, row), GetColumn(xs, col), wrt, block);

                // if we aren't on the diagonal, copy the block of the covariance matrix we just added to the upper triangular part of the covariance matrix
		if(col!=row)
		{
		    for(int j=0; j<coDim; ++j)
		    {
			for(int i=0; i<coDim; ++i)
			{
			    derivs(row*coDim + i, col*coDim + j) = derivs(col*coDim + j, row*coDim + i);
			}
		    }
		}
    	    }
    	}
    }


    
    template<typename PosMatrixType, typename CovMatrixType>
    void FillCovarianceImpl(PosMatrixType const& xs, CovMatrixType & cov) const
    {	
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

    	assert(GetShape(cov, 0) == coDim*GetShape(xs,1));
    	assert(GetShape(cov, 1) == coDim*GetShape(xs,1));

	for(unsigned col=0; col<N1; ++col)
	{
	    for(unsigned row=col; row<N1; ++row)
	    {
		auto block = GetBlock(cov, row*coDim, col*coDim, coDim, coDim);
    		static_cast<const ChildType*>(this)->EvaluateImpl( GetColumn(xs,row), GetColumn(xs, col) , block);

		// if we aren't on the diagonal, copy the block of the covariance matrix we just added to the upper triangular part of the covariance matrix
		if(col!=row)
		{
		    for(int j=0; j<coDim; ++j)
		    {
			for(int i=0; i<coDim; ++i)
			{
			    cov(col*coDim + j, row*coDim + i) = cov(row*coDim + i, col*coDim + j);
			}
		    }
		}
    	    }
    	}
    };

    
    template<typename PosMatrixType, typename CovMatrixType>
    void FillCovarianceImpl(PosMatrixType const& xs, PosMatrixType const& ys, CovMatrixType & cov) const
    {	
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);
    	unsigned N2    = GetShape(ys, 1);

    	assert(GetShape(xs, 0)  == GetShape(ys,0));
    	assert(GetShape(cov, 0) == coDim*GetShape(xs,1));
    	assert(GetShape(cov, 1) == coDim*GetShape(ys,1));

	for(unsigned col=0; col<N2; ++col)
	{
	    for(unsigned row=0; row<N1; ++row)
	    {
		auto block = GetBlock(cov, row*coDim, col*coDim, coDim, coDim);
    		static_cast<const ChildType*>(this)->EvaluateImpl( GetColumn(xs,row), GetColumn(ys,col), block);
    	    }
    	}
	
    };

};


}
}

#endif
