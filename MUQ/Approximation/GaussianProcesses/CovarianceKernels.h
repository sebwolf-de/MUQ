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
template<class ChildType>
class KernelBase
{

public:
    template<typename RightType>
    ProductKernel<ChildType, RightType> operator*(RightType const& kernel2) const
    {
	return ProductKernel<ChildType,RightType>(*reinterpret_cast<const ChildType*>(this),kernel2);
    }

    ChildType operator*(double scale) const
    {
	ChildType newKernel(*reinterpret_cast<const ChildType*>(this));
	newKernel.sigma2 *= scale;
	return newKernel;
    }

    template<typename RightType>
    SumKernel<ChildType,RightType> operator+(RightType const& kernel2) const
    {
	return SumKernel<ChildType,RightType>(*reinterpret_cast<const ChildType*>(this),kernel2);
    }

    template<typename PosMatrixType, typename DerivMatrixType>
    void BuildDerivativeMatrix(PosMatrixType const& xs, int wrt, DerivMatrixType &derivs)
    {
	unsigned codim = ChildType::GetCodim();
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

    	assert(ChildType::GetDim() == GetShape(xs,0));
    	assert(GetShape(derivs, 0) == codim*GetShape(xs,1));
    	assert(GetShape(derivs, 1) == codim*GetShape(xs,1));

	for(unsigned col=0; col<N1; ++col)
	{
	    for(unsigned row=col; row<N1; ++row)
	    {
    		derivs(row, col) = static_cast<ChildType*>(this) -> GetDerivative( GetColumn(xs,row), GetColumn(xs, col), wrt );
		derivs(col, row) = derivs(row,col);
    	    }
    	}
    }

    
    template<typename PosMatrixType, typename CovMatrixType>
    void BuildCovariance(PosMatrixType const& xs, CovMatrixType & cov)
    {
        unsigned codim = ChildType::GetCodim();
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

    	assert(ChildType::GetDim() == GetShape(xs,0));
    	assert(GetShape(cov, 0) == codim*GetShape(xs,1));
    	assert(GetShape(cov, 1) == codim*GetShape(xs,1));

	for(unsigned col=0; col<N1; ++col)
	{
	    for(unsigned row=col; row<N1; ++row)
	    {
    		cov(row, col) = (*static_cast<ChildType*>(this))( GetColumn(xs,row), GetColumn(xs, col) );
		cov(col, row) = cov(row,col);
    	    }
    	}
    };
    
    template<typename PosMatrixType, typename CovMatrixType>
    void BuildCovariance(PosMatrixType const& xs, PosMatrixType const& ys, CovMatrixType & cov)
    {
        unsigned codim = ChildType::GetCodim();
    	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);
    	unsigned N2    = GetShape(ys, 1);

    	assert(ChildType::GetDim() == GetShape(xs,0));
    	assert(GetShape(xs, 0)  == GetShape(ys,0));
    	assert(GetShape(cov, 0) == codim*GetShape(xs,1));
    	assert(GetShape(cov, 1) == codim*GetShape(ys,1));

	for(unsigned col=0; col<N2; ++col)
	{
	    for(unsigned row=0; row<N1; ++row)
	    {
    		cov(row, col) = (*static_cast<ChildType*>(this))( GetColumn(xs,row), GetColumn(ys,col) );
    	    }
    	}
	
    };


    Eigen::MatrixXd GetParamBounds() const
    {
	return paramBounds;
    };
	
protected:

    Eigen::MatrixXd paramBounds;
};


/**
@class TensorProductKernel

Create a covariance kernel with multiple inputs by multiplying several univariate kernels together.
\f[
k([x_1,x_2,\ldots, x_d], [y_1,y_2, \ldots, y_d]) = k_1(x_1,y_1) * k_2(x_2,y_2) * \ldots * k_d(x_d,y_d)
\f]

 */
template<class... KTypes>
class TensorProductKernel : public KernelBase<TensorProductKernel<KTypes...>>
{
public:
    
    TensorProductKernel(const KTypes&... kernelsIn) : kernels(std::make_tuple(kernelsIn...)){};

    template<typename PosMatrixType1, typename PosMatrixType2>
    double operator()(PosMatrixType1 const& xs, PosMatrixType2 const& ys)
    {
	const unsigned dim = GetShape(xs, 0);
	
	KernelEvaluator<0,std::tuple_size<std::tuple<KTypes...>>::value-1,KTypes...> evaluator;
	Eigen::VectorXd kernelVals(dim);
	
	evaluator.TensorEvaluate(kernels, xs, ys, 0, 0, kernelVals);
	return kernelVals.prod();
    }

    template<typename PosMatrixType1, typename PosMatrixType2>
    double GetDerivative(PosMatrixType1 const& xs, PosMatrixType2 const& ys, int wrt) const
    {
	assert(false);
	return 0.0;
    }

    static unsigned GetCodim()
    {
	return 1;
    };
    
    static unsigned GetDim()
    {
	return std::tuple_size<std::tuple<KTypes...>>::value;
    };

    static unsigned GetNumParams()
    {
	return ParamsGrabber<KTypes...>::GetNumParams();
    };


    static unsigned GetNumConstraints()
    {
	return GetDim()-1 + ParamsGrabber<KTypes...>::GetNumConstraints();
    };

    Eigen::VectorXd GetParams() const
    {
	assert(false);
	return Eigen::VectorXd();
    }

    
    template<typename Derived>   
    void SetParams(Eigen::DenseCoeffsBase<Derived> const& params)
    {
	KernelEvaluator<0,std::tuple_size<std::tuple<KTypes...>>::value-1,KTypes...>::SetParams(kernels,params);
    }

	
private:
    
    std::tuple<KTypes...> kernels;
};


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
    
    WhiteNoiseKernel(const double sigma2In,
		     const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0) = sigmaBounds[0]; // lower bound on sigma2
	paramBounds(1) = sigmaBounds[1]; // upper bound on sigma2
    };

    template<typename VecType1, typename VecType2>
    double operator()(VecType1 const& x1, VecType2 const& x2 ) const
    {
	return (*this)(x1(0),x2(0));
    }

    template<typename VecType1, typename VecType2>
    double GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt) const
    {
	assert(wrt==0);

	return (std::abs<double>(x1(0)-x2(0))<1e-14) ? 1.0 : 0.0;
    }
	
    double operator()(const double x1, const double x2) const
    {
	return (std::abs<double>(x1-x2)<1e-14) ? sigma2 : 0.0;
    }

    static unsigned GetNumParams()
    {
	return 1;
    };

    static unsigned GetNumConstraints()
    {
	return 0;
    }


    
    Eigen::VectorXd GetParams() const
    {
	return sigma2*Eigen::VectorXd::Ones(1);
    }
    
    template<typename VecType>
    void SetParams(VecType const& params)
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
    SquaredExpKernel(const double sigma2In,
		     const double lengthIn,
	             const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	             const Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()}) : sigma2(sigma2In), length(lengthIn)
    {
	paramBounds.resize(2,2);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
	
	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);
    };

    template<typename VecType1, typename VecType2>
    double operator()(VecType1 const& x1, VecType2 const& x2 ) const
    {
	return (*this)(x1(0),x2(0));
    }

    template<typename VecType1, typename VecType2>
    double GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt) const
    {
	assert(wrt<GetNumParams());

	if(wrt==0) // derivative wrt sigma2
	{
	    return exp(-0.5*pow((x1(0)-x2(0))/length,2.0));
	}
	else if(wrt==1) // derivative wrt length
	{
	    return sigma2 * exp(-0.5*pow((x1(0)-x2(0))/length,2.0)) * pow(x1(0)-x2(0),2.0) * pow(length, -3.0);
	}
	else
	{
	    assert(false);
	    return 0.0;
	}
    }

    double operator()(const double x1, const double x2) const
    {
	return sigma2*exp(-0.5*pow((x1-x2)/length,2.0));
    }

    static unsigned GetDim()
    {
	return 1;
    }
    
    static unsigned GetCodim()
    {
	return 1;
    };
    
    static unsigned GetNumParams()
    {
	return 2;
    };

    static unsigned GetNumConstraints()
    {
	return 0;
    }
	
    Eigen::VectorXd GetParams() const
    {
	return Eigen::Vector2d{sigma2,length};
    }
    
    template<typename VecType>
    void SetParams(VecType const& params)
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
    ConstantKernel(const double sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : sigma2(sigma2In)
    {
	paramBounds.resize(2,1);
	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);
    };

    template<typename VecType>
    double operator()(VecType const& x1, VecType const& x2 ) const
    {
	return (*this)(x1(0),x2(0));
    }

    template<typename VecType1, typename VecType2>
    double GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt) const
    {
	assert(wrt==0);

	return 1.0;
    }


    double operator()(const double x1, const double x2) const
    {
	return sigma2;
    }

    static unsigned GetDim()
    {
	return 1;
    }
    
    static unsigned GetCodim()
    {
	return 1;
    };
    
    static unsigned GetNumParams()
    {
	return 1;
    };

    static unsigned GetNumConstraints()
    {
	return 0;
    }
	
    Eigen::VectorXd GetParams() const
    {
	return sigma2*Eigen::VectorXd::Ones(1);
    }

    template<typename VecType>
    void SetParams(VecType const& params)
    {
        sigma2 = params(0);
    }

private:
    double sigma2;
};

/** 

@class CoregionalKernel
This kernel supports coregionalization for modeling vector-valued Gaussian processes.

 */
template<class... KTypes>
class CoregionalKernel : public KernelBase<CoregionalKernel<KTypes...>>
{

public:
    CoregionalKernel(Eigen::MatrixXd const& coCov, const KTypes&... kernelsIn) : kernels(std::make_tuple(kernelsIn...))
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(coCov);
        A = eigSolver.eigenvectors()*eigSolver.eigenvalues().cwiseSqrt().asDiagonal();
    };

    template<typename PosMatrixType, typename CovMatrixType>
    void BuildCovariance(PosMatrixType const& xs, PosMatrixType const& ys, CovMatrixType & cov)
    {
        unsigned codim = GetCodim();
	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);
	unsigned N2    = GetShape(ys, 1);

	assert(GetDim() == GetShape(xs,0));
	assert(GetShape(xs, 0)  == GetShape(ys,0));
	assert(GetShape(cov, 0) == codim*GetShape(xs,1));
	assert(GetShape(cov, 1) == codim*GetShape(ys,1));

        Eigen::VectorXd rhoVec(codim);
	for(unsigned row=0; row<N1; ++row)
	{
            for(unsigned col=0; col<N2; ++col)
	    {        
                KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::Evaluate(kernels, xs, ys, row, col, rhoVec);

	        cov.block(row*codim,col*codim,codim,codim) = A * rhoVec.asDiagonal() * A.transpose(); 
	    }
	}
    };

    template<typename PosMatrixType, typename CovMatrixType>
    void BuildCovariance(PosMatrixType const& xs, CovMatrixType & cov)
    {
        unsigned codim = GetCodim();
	unsigned dim   = GetShape(xs, 0);
        unsigned N1    = GetShape(xs, 1);

	assert(GetDim() == GetShape(xs,0));
	assert(GetShape(cov, 0) == codim*GetShape(xs,1));
	assert(GetShape(cov, 1) == codim*GetShape(xs,1));

        Eigen::VectorXd rhoVec(codim);
	for(unsigned row=0; row<N1; ++row)
	{
            for(unsigned col=0; col<=row; ++col)
	    {        
                KernelEvaluator<0, std::tuple_size<std::tuple<KTypes...>>::value, KTypes...>::Evaluate(kernels, xs, xs, row, col, rhoVec);

	        cov.block(row*codim, col*codim, codim, codim) = A * rhoVec.asDiagonal() * A.transpose();
		cov.block(col*codim, row*codim, codim, codim) = cov.block(row*codim, col*codim, codim, codim).transpose().eval(); 
	    }
	}
    };
	
    static unsigned GetDim()
    {
        return std::tuple_element<0, std::tuple<KTypes...> >::type::GetDim();
    };

    static unsigned GetCodim()
    {
	return std::tuple_size<std::tuple<KTypes...>>::value;;
    };
    
    static unsigned GetNumParams()
    {
        unsigned codim = std::tuple_size<std::tuple<KTypes...>>::value;
        return (codim+1)*(codim)/2 + ParamsGrabber<KTypes...>::GetNumParams();
    };

    static unsigned GetNumConstraints()
    {
	return ParamsGrabber<KTypes...>::GetNumConstraints();
    };
	
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
    PeriodicKernel(const double sigma2In,
		   const double lengthIn,
		   const double periodIn,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()},
	           const Eigen::Vector2d lengthBounds = {1e-10, std::numeric_limits<double>::infinity()},
	           const Eigen::Vector2d periodBounds = {1e-10, std::numeric_limits<double>::infinity()}) : sigma2(sigma2In), length(lengthIn), period(periodIn)
    {
	paramBounds.resize(2,3);

	paramBounds(0,0) = sigmaBounds(0);
	paramBounds(1,0) = sigmaBounds(1);

	paramBounds(0,1) = lengthBounds(0);
	paramBounds(1,1) = lengthBounds(1);

	paramBounds(0,2) = periodBounds(0);
	paramBounds(1,2) = periodBounds(1);
    };

    template<typename VecType>
    double operator()(VecType const& x1, VecType const& x2 ) const
    {
	return (*this)(x1(0),x2(0));
    }


    template<typename VecType1, typename VecType2>
    double GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt) const
    {
	assert(wrt < GetNumParams());

	if( wrt == 0)
	{
	    return exp(-2.0 * pow(sin(pi*(x1(0)-x2(0))/period),2.0) / pow(length,2.0));

	}
	else if(wrt==1)
	{
	    return 4.0 * sigma2 * exp(-2.0 * pow(sin(pi*(x1(0)-x2(0))/period),2.0) / pow(length,2.0)) *  std::pow(sin(pi*(x1(0)-x2(0))/period),2.0) * std::pow(length,-3.0);
	}
	else if(wrt==2)
	{
	    double temp = pi*(x1(0)-x2(0))/period;
	    return 4.0 * sigma2 * exp(-2.0 * pow(sin(temp),2.0) / pow(length,2.0)) * sin(temp) * cos(temp) * pi * (x1(0)-x2(0)) * std::pow(length*period,-2.0);
	}
	else
	{
	    assert(false);
	    return 0.0;
	}
    }

    
    double operator()(const double x1, const double x2) const
    {
	return sigma2 * exp(-2.0 * pow(sin(pi*(x1-x2)/period),2.0) / pow(length,2.0));
    }

    

    static unsigned GetNumParams()
    {
	return 3;
    };

    static unsigned GetNumConstraints()
    {
	return 0;
    }
	
    Eigen::VectorXd GetParams() const
    {
	return Eigen::Vector3d{sigma2, length, period};
    };
    
    template<typename VecType>
    void SetParams(VecType const& params)
    {
        sigma2 = params(0);
	length = params(1);
	period = params(2);
    }

    static unsigned GetDim()
    {
	return 1;
    };
    
    static unsigned GetCodim()
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
    ProductKernel(LeftType const& kernel1In, RightType const& kernel2In) : kernel1(kernel1In), kernel2(kernel2In){};

    template<typename VecType>
    double operator()(VecType const& x1, VecType const& x2 ) const
    {
	return kernel1(x1,x2)*kernel2(x1,x2);
    }

    template<typename VecType1, typename VecType2>
    double GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt) const
    {
	assert(wrt < GetNumParams());

	if(wrt < LeftType::GetNumParams() )
	{

	    return kernel1.GetDerivative(x1, x2, wrt) * kernel2(x1,x2);
	}
	else
	{
	    return kernel1(x1,x2) * kernel2.GetDerivative(x1,x2, wrt-LeftType::GetNumParams());
	}
    }
    
    static unsigned GetDim()
    {
        return LeftType::GetDim();
    };
    
    static unsigned GetNumParams()
    {
	return LeftType::GetNumParams() + RightType::GetNumParams();
    };

    static unsigned GetNumConstraints()
    {
	return 1 + LeftType::GetNumConstraints() + RightType::GetNumConstraints();
    }

    Eigen::MatrixXd GetParamBounds() const
    {
	Eigen::MatrixXd bounds(2,GetNumParams());

	bounds.block(0,0,2,kernel1.GetNumParams()) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.GetNumParams(),2,kernel2.GetNumParams()) = kernel2.GetParamBounds();

	return bounds;
    };
	
    Eigen::VectorXd GetParams() const
    {
	Eigen::VectorXd params(GetNumParams());

	params.head(kernel1.GetNumParams()) = kernel1.GetParams();
	params.tail(kernel2.GetNumParams()) = kernel2.GetParams();

	return params;
    };
    
    void SetParams(Eigen::VectorXd const& params)
    {
        kernel1.SetParams(params.head(kernel1.GetNumParams()));
	kernel2.SetParams(params.tail(kernel2.GetNumParams()));
    }

    static unsigned GetCodim()
    {
	return LeftType::GetCodim();
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
    SumKernel(LeftType const& kernel1In, RightType const& kernel2In) : kernel1(kernel1In), kernel2(kernel2In){};

    template<typename VecType>
    double operator()(VecType const& x1, VecType  const& x2 ) const
    {
        return kernel1(x1,x2) + kernel2(x1,x2);
    }

    template<typename VecType1, typename VecType2>
    double GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt) const
    {
	assert(wrt < GetNumParams());

	if(wrt < LeftType::GetNumParams() )
	{

	    return kernel1.GetDerivative(x1, x2, wrt);
	}
	else
	{
	    return kernel2.GetDerivative(x1,x2, wrt-LeftType::GetNumParams());
	}
    }
    
    static unsigned GetDim()
    {
        return LeftType::GetDim();
    };

    static unsigned GetCodim()
    {
        return LeftType::GetCodim();
    };
    
    static unsigned GetNumParams()
    {
	return LeftType::GetNumParams() + RightType::GetNumParams();
    };

    static unsigned GetNumConstraints()
    {
	return LeftType::GetNumConstraints() + RightType::GetNumConstraints();
    }

    Eigen::MatrixXd GetParamBounds() const
    {
	Eigen::MatrixXd bounds(2,GetNumParams());

	bounds.block(0,0,2,kernel1.GetNumParams()) = kernel1.GetParamBounds();
	bounds.block(0,kernel1.GetNumParams(),2,kernel2.GetNumParams()) = kernel2.GetParamBounds();

	return bounds;
    };
	
    Eigen::VectorXd GetParams() const
    {
	Eigen::VectorXd params(GetNumParams());

	params.head(kernel1.GetNumParams()) = kernel1.GetParams();
	params.tail(kernel2.GetNumParams()) = kernel2.GetParams();

	return params;
    };
        
    void SetParams(Eigen::VectorXd const& params)
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
template<class... KTypes>
TensorProductKernel<KTypes...> TensorTie(const KTypes&... kernels)
{
    return TensorProductKernel<KTypes...>(kernels...);
};


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
