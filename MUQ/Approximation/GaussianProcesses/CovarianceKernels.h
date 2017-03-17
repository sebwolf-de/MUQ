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
class KernelBase
{

    static std::vector<unsigned> BuildDimInds(unsigned dim)
    {
	std::vector<unsigned> output(dim);
	for(int i=0; i<dim; ++i)
	    output[i] = i;

	return output;
    }
    
public:

    KernelBase(unsigned dim) : KernelBase(BuildDimInds(dim)){};
    KernelBase(std::vector<unsigned> dimIndsIn) : dimInds(dimIndsIn){};

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
    
    
    /* KernelBase operator*(double scale) const */
    /* { */
    /* 	ChildType newKernel(*reinterpret_cast<const ChildType*>(this)); */
    /* 	newKernel.sigma2 *= scale; */
    /* 	return newKernel; */
    /* } */

    template<typename PosMatrixType, typename DerivMatrixType>
    void BuildDerivativeMatrix(PosMatrixType const& xs, int wrt, DerivMatrixType &derivs)
    {
	unsigned codim = GetCodim();
	assert(codim==1);
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

    	assert(GetShape(derivs, 0) == codim*GetShape(xs,1));
    	assert(GetShape(derivs, 1) == codim*GetShape(xs,1));

	for(unsigned col=0; col<N1; ++col)
	{
	    for(unsigned row=col; row<N1; ++row)
	    {
		auto block = GetBlock(derivs, row*codim, col*codim, codim, codim);
    		static_cast<ChildType*>(this)->GetDerivative( GetColumn(xs, row), GetColumn(xs, col), wrt, block);

                // if we aren't on the diagonal, copy the block of the covariance matrix we just added to the upper triangular part of the covariance matrix
		if(col!=row)
		{
		    for(int j=0; j<codim; ++j)
		    {
			for(int i=0; i<codim; ++i)
			{
			    derivs(row*codim + i, col*codim + j) = derivs(col*codim + j, row*codim + i);
			}
		    }
		}
    	    }
    	}
    }

    
    template<typename PosMatrixType, typename CovMatrixType>
    void BuildCovariance(PosMatrixType const& xs, CovMatrixType & cov)
    {
        unsigned codim = GetCodim();
	assert(codim==1);
	
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

    	assert(GetShape(cov, 0) == codim*GetShape(xs,1));
    	assert(GetShape(cov, 1) == codim*GetShape(xs,1));

	for(unsigned col=0; col<N1; ++col)
	{
	    for(unsigned row=col; row<N1; ++row)
	    {
		auto block = GetBlock(cov, row*codim, col*codim, codim, codim);
    		static_cast<ChildType*>(this)->Evaluate( GetColumn(xs,row), GetColumn(xs, col) , block);

		// if we aren't on the diagonal, copy the block of the covariance matrix we just added to the upper triangular part of the covariance matrix
		if(col!=row)
		{
		    for(int j=0; j<codim; ++j)
		    {
			for(int i=0; i<codim; ++i)
			{
			    cov(row*codim + i, col*codim + j) = cov(col*codim + j, row*codim + i);
			}
		    }
		}
    	    }
    	}
    };
    
    template<typename PosMatrixType, typename CovMatrixType>
    void BuildCovariance(PosMatrixType const& xs, PosMatrixType const& ys, CovMatrixType & cov)
    {
        unsigned codim = GetCodim();
	assert(codim==1);
	
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);
    	unsigned N2    = GetShape(ys, 1);

    	assert(GetShape(xs, 0)  == GetShape(ys,0));
    	assert(GetShape(cov, 0) == codim*GetShape(xs,1));
    	assert(GetShape(cov, 1) == codim*GetShape(ys,1));

	for(unsigned col=0; col<N2; ++col)
	{
	    for(unsigned row=0; row<N1; ++row)
	    {
		auto block = GetBlock(cov, row*codim, col*codim, codim, codim);
    		static_cast<ChildType*>(this)->Evaluate( GetColumn(xs,row), GetColumn(ys,col), block);
    	    }
    	}
	
    };


    virtual Eigen::MatrixXd GetParamBounds() const
    {
	return paramBounds;
    };


    virtual int GetDim() const{return dimInds.size();};
    virtual int GetCodim() const{return 1;};
    virtual int GetNumParams() const{return 0;};
    virtual Eigen::VectorXd GetParams() const{return Eigen::VectorXd();};

    virtual void SetParams(Eigen::VectorXd const& params){};


    const std::vector<unsigned> dimInds;

protected:

    
    Eigen::MatrixXd paramBounds;
};


/* /\** */
/* @class TensorProductKernel */

/* Create a covariance kernel with multiple inputs by multiplying several univariate kernels together. */
/* \f[ */
/* k([x_1,x_2,\ldots, x_d], [y_1,y_2, \ldots, y_d]) = k_1(x_1,y_1) * k_2(x_2,y_2) * \ldots * k_d(x_d,y_d) */
/* \f] */

/*  *\/ */
/* template<class... KTypes> */
/* class TensorProductKernel : public KernelBase<TensorProductKernel<KTypes...>> */
/* { */
/* public: */
    
/*     TensorProductKernel(const KTypes&... kernelsIn) : kernels(std::make_tuple(kernelsIn...)) */
/*     { */
/* 	dim = std::tuple_size<std::tuple<KTypes...>>::value; */
/*     }; */

/*     template<typename PosMatrixType1, typename PosMatrixType2> */
/*     double operator()(PosMatrixType1 const& xs, PosMatrixType2 const& ys) */
/*     { */
/* 	const unsigned dim = GetShape(xs, 0); */
	
/* 	KernelEvaluator<0,std::tuple_size<std::tuple<KTypes...>>::value-1,KTypes...> evaluator; */
/* 	Eigen::VectorXd kernelVals(dim); */
	
/* 	evaluator.TensorEvaluate(kernels, xs, ys, 0, 0, kernelVals); */
/* 	return kernelVals.prod(); */
/*     } */

/*     template<typename PosMatrixType1, typename PosMatrixType2> */
/*     double GetDerivative(PosMatrixType1 const& xs, PosMatrixType2 const& ys, int wrt) const */
/*     { */
/* 	assert(false); */
/* 	return 0.0; */
/*     } */

/*     virtual int GetCodim() const override */
/*     { */
/* 	return 1; */
/*     }; */
    
/*     virtual int GetNumParams() const override */
/*     { */
/* 	return ParamsGrabber<KTypes...>::GetNumParams(); */
/*     }; */

/*     virtual Eigen::VectorXd GetParams() const override */
/*     { */
/* 	assert(false); */
/* 	return Eigen::VectorXd(); */
/*     } */

/*     virtual void SetParams(Eigen::VectorXd const& params) override */
/*     { */
/* 	KernelEvaluator<0,std::tuple_size<std::tuple<KTypes...>>::value-1,KTypes...>::SetParams(kernels,params); */
/*     } */

	
/* private: */
    
/*     std::tuple<KTypes...> kernels; */
/* }; */


/**

@class WhiteNoiseKernel

This class implements a kernel of the form
\f[
k(x,y) = \var^2\delta(x,y)
\f]
where \f$\delta(x,y)=1\f$ if \f$x=y\f$ and \f$0\f$ otherwise.

 */
class WhiteNoiseKernel : public KernelBase<WhiteNoiseKernel>
{

public:

    template<typename dimType>
    WhiteNoiseKernel(dimType dim,
		     const double sigma2In,
		     const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelBase<WhiteNoiseKernel>(dim), sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0) = sigmaBounds[0]; // lower bound on sigma2
	paramBounds(1) = sigmaBounds[1]; // upper bound on sigma2
    };

    template<typename VecType1, typename VecType2, typename MatrixType>
    void Evaluate(VecType1 const& x1, VecType2 const& x2, MatrixType &cov) const
    {
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	cov(0,0) = (dist<1e-14) ? sigma2 : 0.0;
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType &derivs) const
    {
	assert(wrt==0);
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));

	derivs(0,0) = (dist<1e-14) ? 1.0 : 0.0;
    }
	

    virtual int GetNumParams() const override
    {
	return 1;
    };

    
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
class SquaredExpKernel : public KernelBase<SquaredExpKernel>
{

public:

    template<typename dimType>
    SquaredExpKernel(dimType dimIn,
		     const double sigma2In,
		     const double lengthIn,
	             const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	             const Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelBase<SquaredExpKernel>(dimIn), sigma2(sigma2In), length(lengthIn)
    {
	paramBounds.resize(2,2);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
	
	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);
    };

    template<typename VecType1, typename VecType2, typename MatrixType>
    void Evaluate(VecType1 const& x1, VecType2 const& x2, MatrixType &cov ) const
    {
	assert(GetShape(x2,0)==GetShape(x1,0));

	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	
	cov(0,0) = sigma2*exp(-0.5*pow(dist/length,2.0));
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	assert(wrt<GetNumParams());

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
    
    virtual int GetCodim() const override
    {
	return 1;
    };
    
    virtual int GetNumParams() const override
    {
	return 2;
    };
	
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
class ConstantKernel : public KernelBase<ConstantKernel>
{

public:

    template<typename dimType>
    ConstantKernel(dimType dim,
	           const double sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : ConstantKernel(dim, sigma2In*Eigen::MatrixXd::Ones(1,1), sigmaBounds){};
	
    template<typename dimType>
    ConstantKernel(dimType dim,
	           Eigen::MatrixXd const& sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelBase<ConstantKernel>(dim), sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
    };

    template<typename VecType, typename MatrixType>
    void Evaluate(VecType const& x1, VecType const& x2, MatrixType & cov ) const
    {
	cov = sigma2;
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
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

    // The parameters in this kernel are the lower triangular components of the covariance.
    virtual int GetNumParams() const override
    {
	return sigma2.rows()*(sigma2.cols()+1)/2;
    };
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd output(GetNumParams());
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
};

/** 

@class CoregionalKernel
This kernel supports coregionalization for modeling vector-valued Gaussian processes.

 */
template<class... KTypes>
class CoregionalKernel : public KernelBase<CoregionalKernel<KTypes...>>
{

public:
    
    CoregionalKernel(Eigen::MatrixXd const& Gamma, const KTypes&... kernelsIn) : KernelBase<CoregionalKernel<KTypes...>>(Gamma.rows()), kernels(std::make_tuple(kernelsIn...))
    {
	// Make sure the matrix and kernels are the same size
	assert(Gamma.cols()==std::tuple_size<std::tuple<KTypes...>>::value);

	// Compute the eigenvalue decomposition of the covariance Gamma to get the matrix A
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver;
	eigSolver.compute(Gamma, true);
	A = eigSolver.eigenvectors()*eigSolver.eigenvalues().cwiseSqrt().asDiagonal();
    };

    template<typename VecType, typename MatrixType>
    void Evaluate(VecType const& x1, VecType const& x2, MatrixType & cov ) const
    {
	Eigen::VectorXd rhoVec(GetCodim());
	KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::Evaluate(kernels, x1, x2, rhoVec);
	cov = A * rhoVec.asDiagonal() * A.transpose();
    }

    template<typename VecType1, typename VecType2, typename MatrixType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
    {
	int kernelInd = 0;
	double kernelDeriv = KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::GetDeriv(kernels, x1, x2, wrt, 0, kernelInd);
	derivs = Eigen::MatrixXd(A.col(kernelInd) * kernelDeriv * A.col(kernelInd).transpose());
    }

    
    virtual int GetCodim() const override
    {
	return std::tuple_size<std::tuple<KTypes...>>::value;;
    };
    
    virtual int GetNumParams() const override
    {
        unsigned codim = std::tuple_size<std::tuple<KTypes...>>::value;
        return ParamsGrabber<KTypes...>::GetNumParams();
    };

    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(GetNumParams());
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

};





/** 

@class PeriodicKernel

\f[
k(x,y) = \exp\left[-\frac{2}{L^2} \sin^2\left(\frac{\pi}{p}(x-y)\right)\right]

 */
class PeriodicKernel : public KernelBase<PeriodicKernel>
{

public:

    template<typename dimType>
    PeriodicKernel(dimType dim,
		   const double sigma2In,
		   const double lengthIn,
		   const double periodIn,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	           const Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()},
	           const Eigen::Vector2d periodBounds = {1e-10, std::numeric_limits<double>::infinity()}) : KernelBase<PeriodicKernel>(dim), sigma2(sigma2In), length(lengthIn), period(periodIn)
    {
	paramBounds.resize(2,3);

	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);

	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);

	paramBounds(0,2) = periodBounds(0);
	paramBounds(1,2) = periodBounds(1);
    };

    template<typename VecType, typename MatType>
    void Evaluate(VecType const& x1, VecType const& x2, MatType & cov) const
    {
	double dist = CalcDistance(GetSlice(x1, dimInds), GetSlice(x2, dimInds));
	cov(0,0) = sigma2 * exp(-2.0 * pow(sin(pi*dist/period),2.0) / pow(length,2.0));
    }


    template<typename VecType1, typename VecType2, typename MatType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < GetNumParams());

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

    virtual int GetNumParams() const override
    {
	return 3;
    };

	
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

    
    virtual int GetCodim() const override
    {
	return 1;
    };
    
private:
    double sigma2;
    double length;
    double period;

    const double pi = 4.0 * atan(1.0); //boost::math::constants::pi<double>();
};


/**

@class ProductKernel

\f[
k(x,y) = k_1(x,y)*k_2(x,y)
\f]

 */
template<typename LeftType, typename RightType>
class ProductKernel : public KernelBase<ProductKernel<LeftType,RightType>>
{

public:
    ProductKernel(LeftType const& kernel1In, RightType const& kernel2In) : KernelBase<ProductKernel<LeftType,RightType>>(kernel1In.GetDim()), kernel1(kernel1In), kernel2(kernel2In)
    {
	assert(kernel1.GetCodim()==kernel2.GetCodim());
    };

    template<typename VecType, typename MatType>
    void Evaluate(VecType const& x1, VecType const& x2, MatType & cov ) const
    {
	Eigen::MatrixXd temp1(GetCodim(), GetCodim());
	Eigen::MatrixXd temp2(GetCodim(), GetCodim());

	kernel1.Evaluate(x1,x2, temp1);
	kernel2.Evaluate(x1,x2, temp2);

	cov = Eigen::MatrixXd(temp1.array() * temp2.array());
    }

    template<typename VecType1, typename VecType2, typename MatType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < GetNumParams());

	Eigen::MatrixXd temp1(GetCodim(), GetCodim());
	Eigen::MatrixXd temp2(GetCodim(), GetCodim());

	if(wrt < kernel1.GetNumParams() )
	{	    
	    kernel1.GetDerivative(x1, x2, wrt, temp1);
	    kernel2.Evaluate(x1,x2, temp2);
	}
	else
	{
	    kernel1.Evaluate(x1,x2, temp1);
	    kernel2.GetDerivative(x1,x2, wrt-kernel1.GetNumParams(), temp2);
	}
	derivs = Eigen::MatrixXd(temp1.array() * temp2.array());
    }
    
    virtual int GetNumParams() const override
    {
	return kernel1.GetNumParams() + kernel2.GetNumParams();
    };

    virtual Eigen::MatrixXd GetParamBounds() const override
    {
	Eigen::MatrixXd bounds(2,GetNumParams());

	bounds.block(0,0,2,kernel1.GetNumParams()) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.GetNumParams(),2,kernel2.GetNumParams()) = kernel2.GetParamBounds();

	return bounds;
    };
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(GetNumParams());

	params.head(kernel1.GetNumParams()) = kernel1.GetParams();
	params.tail(kernel2.GetNumParams()) = kernel2.GetParams();

	return params;
    };
    
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        kernel1.SetParams(params.head(kernel1.GetNumParams()));
	kernel2.SetParams(params.tail(kernel2.GetNumParams()));
    }

    virtual int GetCodim() const override
    {
	return kernel1.GetCodim();
    };
    
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
class SumKernel : public KernelBase<SumKernel<LeftType,RightType>>
{

public:
    SumKernel(LeftType const& kernel1In, RightType const& kernel2In) : KernelBase<SumKernel<LeftType,RightType>>(kernel1In.GetDim()), kernel1(kernel1In), kernel2(kernel2In)
    {
	assert(kernel1.GetCodim() == kernel2.GetCodim());
    };

    template<typename VecType, typename MatType>
	void Evaluate(VecType const& x1, VecType  const& x2 , MatType & cov) const
    {
	Eigen::MatrixXd temp1(GetCodim(), GetCodim());
	Eigen::MatrixXd temp2(GetCodim(), GetCodim());
	
        kernel1.Evaluate(x1,x2, temp1);
	kernel2.Evaluate(x1,x2, temp2);
	cov = temp1 + temp2;
    }

    template<typename VecType1, typename VecType2, typename MatType>
    void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
    {
	assert(wrt < GetNumParams());

        if(wrt < kernel1.GetNumParams() )
	{
	    return kernel1.GetDerivative(x1, x2, wrt, derivs);
	}
	else
	{
	    return kernel2.GetDerivative(x1,x2, wrt-kernel1.GetNumParams(), derivs);
	}
    }
    

    virtual int GetCodim() const override
    {
        return kernel1.GetCodim();
    };
    
    virtual int GetNumParams() const override
    {
	return kernel1.GetNumParams() + kernel2.GetNumParams();
    };


    virtual Eigen::MatrixXd GetParamBounds() const override
    {
	Eigen::MatrixXd bounds(2,GetNumParams());

	bounds.block(0,0,2,kernel1.GetNumParams()) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.GetNumParams(),2,kernel2.GetNumParams()) = kernel2.GetParamBounds();

	return bounds;
    };
	
    virtual Eigen::VectorXd GetParams() const override
    {
	Eigen::VectorXd params(GetNumParams());

	params.head(kernel1.GetNumParams()) = kernel1.GetParams();
	params.tail(kernel2.GetNumParams()) = kernel2.GetParams();

	return params;
    };
        
    virtual void SetParams(Eigen::VectorXd const& params) override
    {
        kernel1.SetParams(params.head(kernel1.GetNumParams()));
	kernel2.SetParams(params.tail(kernel2.GetNumParams()));
    }

    LeftType  kernel1;
    RightType kernel2;
};


/** 

Return the tensor product of several scalar valued kernels.

USAGE:
   TensorTie(k1,k2)
   TensorTie(k1,k2,k3)
   TensorTie(k1,k2,k3,k4)

*/
/* template<class... KTypes> */
/* TensorProductKernel<KTypes...> TensorTie(const KTypes&... kernels) */
/* { */
/*     return TensorProductKernel<KTypes...>(kernels...); */
/* }; */


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
