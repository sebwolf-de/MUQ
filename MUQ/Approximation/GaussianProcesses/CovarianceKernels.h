#ifndef COVARIANCEKERNELS_H_
#define COVARIANCEKERNELS_H_

#include <assert.h>
#include <memory>
#include <iostream>
#include <fstream>

#include "MUQ/Approximation/GaussianProcesses/Utilities.h"

namespace muq
{
namespace Approximation
{

/**
\defgroup Kernels Covariance Kernels

Gaussian processes are defined by their mean function and covariance kernel.  In this group, we provide tools for constructing and using covariance kernels.  Several simple kernels are provided as well as tools for combining simple kernels into more complicated kernels.

C++11 introduced the "auto" keyword that tells the compiler to try to automatically detect the type of a variable.  For example, consider a kernel of the form
\f[
k(x,y) = k_1(x,y) * k_2(x,y) + k_3(x,y),
\f]
where \f$k_1\f$ is a squared exponential kernel, \f$k_2\f$ is a periodic kernel and \f$k_3\f$ is a white noise kernel.  In code, this can be implemented as
\code{.cpp}
auto kernel1 = SquaredExpKernel(var1, length1);
auto kernel2 = PeriodicKernel(var2, length2, period);
auto kernel3 = WhiteNoiseKernel(var3);

auto kernel = kernel1*kernel2 + kernel3;
\endcode
or, more succinctly, as 
\code{.cpp}
auto kernel = SquaredExpKernel(var1, length1) * PeriodicKernel(var2, length2, period) + WhiteNoiseKernel(var3);
\endcode

In either case, the kernel variable is an instance of "SumKernel<ProductKernel<SquaredExpKernel,PeriodicKernel>, WhiteNoiseKernel>".  The "auto" keyword allows us to avoid typing this long type.

@{
 
*/

template<typename LeftType, typename RightType>
class ProductKernel;

template<typename LeftType, typename RightType>
class SumKernel;


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


    virtual Eigen::MatrixXd Evaluate(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2) const = 0;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::MatrixXd             const& ys,
				Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
				Eigen::Ref<Eigen::MatrixXd>        cov) const = 0;

    virtual void FillDerivativeMatrix(Eigen::MatrixXd             const& xs,
				      unsigned                           wrt,
				      Eigen::Ref<Eigen::MatrixXd>        derivs) const = 0;
    
    
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


/**

\class KernelBase
\brief Base class for covariance kernels
\details This class provides common functionality (such as computing Covariance matrices) for all covariance kernels.  It uses the curiously recurring template pattern and requires that child classes implement the following functions
- double operator(VectorType1,VectorType2)
- unsigned GetDim()
- unsigned GetCodim()
- unsigned GetNumParams()
- unsigned GetNumConstraints()
- void SetParams(VectorType)

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
	static_cast<const ChildType*>(this)->EvaluateImpl(x1,x2,output);
	return output;
    }
    
    /* KernelBase operator*(double scale) const */
    /* { */
    /* 	ChildType newKernel(*reinterpret_cast<const ChildType*>(this)); */
    /* 	newKernel.sigma2 *= scale; */
    /* 	return newKernel; */
    /* } */

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
			    cov(row*coDim + i, col*coDim + j) = cov(col*coDim + j, row*coDim + i);
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

/**

@class WhiteNoiseKernel

This class implements a kernel of the form
\f[
k(x,y) = \var^2\delta(x,y)
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

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov) const
    {
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	cov(0,0) = (dist<1e-14) ? sigma2 : 0.0;
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



/**

@class SquaredExpKernel

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

/**

@class ConstantKernel

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

private:
    Eigen::MatrixXd sigma2;

    static unsigned GetNumParams(Eigen::MatrixXd const& cov)
    {
	return 0.5*cov.rows()*(cov.rows()+1);
    }
};

/** 

@class CoregionalKernel
This kernel supports coregionalization for modeling vector-valued Gaussian processes.

 */
template<class... KTypes>
class CoregionalKernel : public KernelImpl<CoregionalKernel<KTypes...>>
{

public:
    
    CoregionalKernel(unsigned               dim,
		     Eigen::MatrixXd const& Gamma,
		     const KTypes&...       kernelsIn) : KernelImpl<CoregionalKernel<KTypes...>>(dim, Gamma.rows(), GetNumParams(kernelsIn...)),
	                                                 kernels(std::make_tuple(kernelsIn...))
    {
	// Make sure the matrix and kernels are the same size
	assert(Gamma.cols()==std::tuple_size<std::tuple<KTypes...>>::value);

	// Compute the eigenvalue decomposition of the covariance Gamma to get the matrix A
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver;
	eigSolver.compute(Gamma, true);
	A = eigSolver.eigenvectors()*eigSolver.eigenvalues().cwiseSqrt().asDiagonal();
    };

    template<typename VecType, typename MatrixType>
    inline void EvaluateImpl(VecType const& x1, VecType const& x2, MatrixType & cov ) const
    {
	Eigen::VectorXd rhoVec(this->coDim);
	KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::Evaluate(kernels, x1, x2, rhoVec);
	cov = A * rhoVec.asDiagonal() * A.transpose();
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	int kernelInd = 0;
	double kernelDeriv = KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::GetDeriv(kernels, x1, x2, wrt, 0, kernelInd);
	derivs = Eigen::MatrixXd(A.col(kernelInd) * kernelDeriv * A.col(kernelInd).transpose());
    }


    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(this->numParams);
	KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::GetParams(kernels, params);
	return params;
    }

    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::SetParams(kernels, params);
    }
	
private:

    Eigen::MatrixXd A;
    std::tuple<KTypes...> kernels;

    template<class Type1, class... OtherTypes>
    static unsigned GetNumParams(Type1 const& kernel1, const OtherTypes&... otherKernels)
    {
        return kernel1.numParams + GetNumParams(otherKernels...);
    };
    template<class Type1>
    static unsigned GetNumParams(Type1 const& kernel1)
    {
        return kernel1.numParams;
    };


};





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


/**

@class ProductKernel

\f[
k(x,y) = k_1(x,y)*k_2(x,y)
\f]

 */
template<typename LeftType, typename RightType>
class ProductKernel : public KernelImpl<ProductKernel<LeftType,RightType>>
{

public:
ProductKernel(LeftType const& kernel1In, RightType const& kernel2In) : KernelImpl<ProductKernel<LeftType,RightType>>(kernel1In.inputDim,
														     std::max(kernel1In.coDim, kernel2In.coDim),
														     kernel1In.numParams + kernel2In.numParams),
	                                                               kernel1(kernel1In),
	                                                               kernel2(kernel2In)
    {
	assert((kernel1.coDim==kernel2.coDim) | (kernel1.coDim==1) | (kernel2.coDim==1));
    };

    template<typename VecType, typename MatType>
    inline void EvaluateImpl(VecType const& x1, VecType const& x2, MatType & cov ) const
    {
	Eigen::MatrixXd temp1(kernel1.coDim, kernel1.coDim);
	Eigen::MatrixXd temp2(kernel2.coDim, kernel2.coDim);

	kernel1.EvaluateImpl(x1,x2, temp1);
	kernel2.EvaluateImpl(x1,x2, temp2);

	if(kernel1.coDim==kernel2.coDim)
	{
	    cov = Eigen::MatrixXd(temp1.array() * temp2.array());
	}
	else if(kernel1.coDim==1)
	{
	    cov = temp1(0,0)*temp2;
	}
	else if(kernel2.coDim==1)
	{
	    cov = temp2(0,0)*temp1;
	}
	else
	{
	    std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
	    assert(false);
	}
	
    }

    template<typename VecType1, typename VecType2, typename MatType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < this->numParams);

	Eigen::MatrixXd temp1(kernel1.coDim, kernel1.coDim);
	Eigen::MatrixXd temp2(kernel2.coDim, kernel2.coDim);

	if(wrt < kernel1.numParams )
	{	    
	    kernel1.GetDerivative(x1, x2, wrt, temp1);
	    kernel2.EvaluateImpl(x1,x2, temp2);
	}
	else
	{
	    kernel1.EvaluateImpl(x1,x2, temp1);
	    kernel2.GetDerivative(x1,x2, wrt-kernel1.numParams, temp2);
	}

	if(kernel1.coDim==kernel2.coDim)
	{
	    derivs = Eigen::MatrixXd(temp1.array() * temp2.array());
	}
	else if(kernel1.coDim==1)
	{
	    derivs = temp1(0,0) * temp2;
	}
	else if(kernel2.coDim==1)
	{
            derivs = temp2(0,0) * temp1;
	}
	else
	{
	    std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
	    assert(false);
	}
	
    }
    
    virtual Eigen::MatrixXd GetParamBounds() const override
    {
	Eigen::MatrixXd bounds(2,this->numParams);

	bounds.block(0,0,2,kernel1.numParams) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.numParams,2,kernel2.numParams) = kernel2.GetParamBounds();

	return bounds;
    };
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(this->numParams);

	params.head(kernel1.numParams) = kernel1.GetParams();
	params.tail(kernel2.numParams) = kernel2.GetParams();

	return params;
    };
    
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        kernel1.SetParams(params.head(kernel1.numParams));
	kernel2.SetParams(params.tail(kernel2.numParams));
    }
    
protected:
    LeftType  kernel1;
    RightType kernel2;

};

/**

@class SumKernel

\f[
k)x,y) = k_1(x,y) + k_2(x,y)
\f]

 */

template<typename LeftType, typename RightType>
class SumKernel : public KernelImpl<SumKernel<LeftType,RightType>>
{

public:
    SumKernel(LeftType const& kernel1In, RightType const& kernel2In) : KernelImpl<SumKernel<LeftType,RightType>>(kernel1In.inputDim,
														 kernel1In.coDim,
														 kernel1In.numParams + kernel2In.numParams),
	                                                               kernel1(kernel1In),
	                                                               kernel2(kernel2In)
    {
	assert(kernel1.coDim == kernel2.coDim);
    };

    template<typename VecType, typename MatType>
    inline void EvaluateImpl(VecType const& x1, VecType  const& x2 , MatType & cov) const
    {
	Eigen::MatrixXd temp1(this->coDim, this->coDim);
	Eigen::MatrixXd temp2(this->coDim, this->coDim);
	
        kernel1.EvaluateImpl(x1,x2, temp1);
	kernel2.EvaluateImpl(x1,x2, temp2);
	cov = temp1 + temp2;
    }

    template<typename VecType1, typename VecType2, typename MatType>
    inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < this->numParams);

        if(wrt < kernel1.numParams )
	{
	    return kernel1.GetDerivative(x1, x2, wrt, derivs);
	}
	else
	{
	    return kernel2.GetDerivative(x1,x2, wrt-kernel1.numParams, derivs);
	}
    }
    
    virtual Eigen::MatrixXd GetParamBounds() const override
    {
	Eigen::MatrixXd bounds(2,this->numParams);

	bounds.block(0,0,2,kernel1.numParams) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.numParams,2,kernel2.numParams) = kernel2.GetParamBounds();

	return bounds;
    };
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(this->numParams);

	params.head(kernel1.numParams) = kernel1.GetParams();
	params.tail(kernel2.numParams) = kernel2.GetParams();

	return params;
    };
        
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        kernel1.SetParams(params.head(kernel1.numParams));
	kernel2.SetParams(params.tail(kernel2.numParams));
    }

    
private:
    LeftType  kernel1;
    RightType kernel2;
};



/**

@class WhiteNoiseKernel

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

    template<typename VecType1, typename VecType2, typename MatrixType>
    inline void EvaluateImpl(VecType1 const& x1, VecType2 const& x2, MatrixType &cov) const
    {
	Eigen::MatrixXd tempCov(K.coDim, K.coDim);
	K.EvaluateImpl(x1,x2,tempCov);
	
	cov = A * tempCov * A.transpose();
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

/**

Create a coregional kernel
@param[in] coCov The joint covariance of the vector-valued process
@param[in] kernels Covariance kernels for the principal components of the covariance.
 
*/
template<class... KTypes>
CoregionalKernel<KTypes...> CoregionTie(Eigen::MatrixXd const& coCov, const KTypes&... kernels)
{
    return CoregionalKernel<KTypes...>(coCov, kernels...);
};


/** @} */ // end of Doxygen covariance kernel group

} // namespace Approximation
} // namespace muq



#endif // #ifndef COVARIANCEKERNELS_H_
