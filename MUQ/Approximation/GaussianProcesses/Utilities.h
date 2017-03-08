
#ifndef UTILITIES_H_
#define UTILITIES_H_

namespace muq
{
namespace Approximation
{

/** 

Templated function to get the shape of an array.  In Eigen parlance, GetShape(mat,0) is the same as mat.rows() and GetShape(mat,1) is the same as mat.cols().

This function is overloaded to handle general types of arrays.  Currently, Kokkos::View<double**> and Eigen matrices are supported.

*/
template<typename MatrixType>
unsigned GetShape(MatrixType const& mat, unsigned dim)
{
    return mat.dimension(dim);
}

template<typename ScalarType, int rows, int cols>
unsigned GetShape(Eigen::Matrix<ScalarType, rows, cols> const& mat, unsigned dim)
{
    assert(dim<2);
    return dim==0 ? mat.rows() : mat.cols();
}


/** 

This class is used like mat.col(colNum) in Eigen, but can be used with more general matrix types (e.g., Kokkos::View<double**>).  Based on a two dimensional matrix (anything implementing operator()(int,int), this class provides a mechanism for accessing elements in a particular column.

This class is usually not constructed directly, but is created using the "ColumnSlice" function, allowing the compiler to detect the correct template type.

 */
template<typename MatType>
class ColumnSlice
{
public:
    ColumnSlice(MatType const& matrixIn, unsigned colIn) : col(colIn), matrix(matrixIn){};

    double const& operator()(unsigned row) const{return matrix(row,col);};
    double const& operator()(unsigned row, unsigned col2) const{assert(col2==0); return matrix(row,col);};
    
    unsigned dimension(unsigned dim) const
    {
	if(dim>0)
	    return 1;
	else
	    return GetShape(matrix,0);
    };
private:
    const unsigned   col;
    MatType        const& matrix;
};

/** 

Grab a particular column of a matrix and return as an instance of the "ColumnSlice" class. 

*/
template<typename MatType>
ColumnSlice<MatType> GetColumn(MatType const& matrix, unsigned col)
{
    return ColumnSlice<MatType>(matrix,col);
}



/**

The tensor product and coregional kernels store a tuple of other kernels with possible different types.  This structure provides a way for looping over the kernels in those tuples and evaluating each one.  This "looping" behavior is achieved through a recursive template that increments the "StartDim" template parameter.

*/
template<unsigned StartDim, unsigned EndDim, class... KTypes>
struct KernelEvaluator
{
    template<typename MatrixType>
    static void Evaluate(std::tuple<KTypes...> const& kernels,
			 MatrixType            const& xs,
			 MatrixType            const& ys,
			 unsigned                     xcol,
			 unsigned                     ycol,
			 Eigen::VectorXd            & output)
    {
	output(StartDim) = std::get<StartDim>(kernels)( GetColumn(xs,xcol), GetColumn(ys,ycol) );
        KernelEvaluator<StartDim+1, EndDim, KTypes...>::Evaluate(kernels, xs, ys, xcol, ycol, output);
	
    }

    template<typename MatrixType>
    static void TensorEvaluate(std::tuple<KTypes...> const& kernels,
			       MatrixType            const& xs,
			       MatrixType            const& ys,
			       unsigned                     xcol,
			       unsigned                     ycol,
			       Eigen::VectorXd            & output)
    {
        output(StartDim) = std::get<StartDim>(kernels)(xs(StartDim,xcol),ys(StartDim,ycol));
        KernelEvaluator<StartDim+1, EndDim, KTypes...>::Evaluate(kernels, xs, ys, xcol, ycol, output);
    }

    template<typename Derived>
    static void SetParams(std::tuple<KTypes...>     const& kernels,
			  Eigen::DenseCoeffsBase<Derived> const& params,
	                  unsigned                         currInd =0)
    {
        auto kernel = std::get<StartDim>(kernels);
	unsigned numParams = kernel.GetNumParams();
	kernel.SetParams(params.segment(currInd,numParams));
	
        KernelEvaluator<StartDim+1, EndDim, KTypes...>::SetParams(kernels, params, currInd+numParams);
    }
};

/**

Structure to end the recursive evaluations of the KernelEvaluator.

 */
template<unsigned StartDim, class... KTypes>
struct KernelEvaluator<StartDim, StartDim, KTypes...>
{
    template<typename MatrixType>
    static void Evaluate(std::tuple<KTypes...> const& kernels,
			 MatrixType            const& xs,
			 MatrixType            const& ys,
			 unsigned                     xcol,
			 unsigned                     ycol,
			 Eigen::VectorXd            & output)
    {
	// do nothing
    }

    template<typename MatrixType>
    static void TensorEvaluate(std::tuple<KTypes...> const& kernels,
			       MatrixType            const& xs,
			       MatrixType            const& ys,
			       unsigned                     xcol,
			       unsigned                     ycol,
			       Eigen::VectorXd            & output)
    {
	// do nothing
    }

    template<typename Derived>
    static void SetParams(std::tuple<KTypes...>     const& kernels,
			  Eigen::DenseCoeffsBase<Derived> const& params,
	                  unsigned                         currInd = 0)
    {
	// do nothing
    }
};


/**

Like "KernelEvaluator", this class helps loop through the tuples used in the TensorProduct and Coregion kernels.  

Here, we provide methods for computing the total number of parameters and number of constraints by adding the parameters and constraints from each kernel in the tuple.

*/
template<class... KTypes>
struct ParamsGrabber
{
    static unsigned GetNumParams()
    {
	assert(false);
    }

    static unsigned GetNumConstraints()
    {
	assert(false);
    }
};

template<class FirstK, class... KTypes>
struct ParamsGrabber<FirstK, KTypes...>
{
    static unsigned GetNumParams()
    {
	return FirstK::GetNumParams() + ParamsGrabber<KTypes...>::GetNumParams();
    }

    static unsigned GetNumConstraints()
    {
	return FirstK::GetNumConstraints() + ParamsGrabber<KTypes...>::GetNumConstraints();
    }
};

template<>
struct ParamsGrabber<>
{
    static unsigned GetNumParams()
    {
	return 0;
    }

    static unsigned GetNumConstraints()
    {
	return 0;
    }
};

} // namespace Approximation

} // namespace muq

#endif 
