
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <Eigen/Dense>
#include <vector>

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

template<typename Derived>
unsigned GetShape(Eigen::Ref<Derived> const& mat, unsigned dim)
{
    assert(dim<2);
    return dim==0 ? mat.rows() : mat.cols();
}

/** @brief Calculates the distance squared between two points defined by vectors v1 and v2. 
    @details Assumes the vectors are the same size and recursively compute the squared distance
             between them.  The recursion is used for numerical accuracy.
*/
template<typename VectorType1, typename VectorType2>
double CalcSquaredDistance(VectorType1 const& v1, VectorType2 const& v2, int startDim=0, int endDim=-1)
{
    if(endDim==-1)
    {
	endDim = GetShape(v1,0);
    }

    const int dim = endDim-startDim;
    const int minDim = 10;
    
    // If the dimension is small enough, just compute the some with a for loop
    if(dim<minDim)
    {
	double output = 0.0;
	for(int i=0; i<dim; ++i)
	{
	    output += std::pow(v1(startDim+i)-v2(startDim+i), 2.0);
	}
	return output;
    }
    else
    {

	int midDim = startDim + std::floor(0.5*dim);
	return CalcSquaredDistance(v1,v2, startDim, midDim) + CalcSquaredDistance(v1,v2, midDim, endDim);
    }

}


/** Calculates the distance between two points defined by vectors v1 and v2. */
template<typename VectorType1, typename VectorType2>
double CalcDistance(VectorType1 const& v1, VectorType2 const& v2)
{
    // Make sure the vectors are the same size
    const int dim = GetShape(v1,0);
    assert(dim==GetShape(v2,0));

    return std::sqrt(CalcSquaredDistance(v1,v2));
}

/** 

This class is used like mat.col(colNum) in Eigen, but can be used with more general matrix types (e.g., Kokkos::View<double**>).  Based on a two dimensional matrix (anything implementing operator()(int,int), this class provides a mechanism for accessing elements in a particular column.

This class is usually not constructed directly, but is created using the "ColumnSlice" function, allowing the compiler to detect the correct template type.

 */
template<typename MatType>
class ColumnSlice
{
public:
    //ColumnSlice(MatType const& matrixIn, unsigned colIn) : col(colIn), matrix(matrixIn){};
    ColumnSlice(MatType const& matrixIn, unsigned colIn) : col(colIn), matrix(matrixIn){};

    double const& operator()(unsigned row) const{return matrix(row, col);};
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

template<typename VecType>
class VectorSlice
{
public:
    VectorSlice(VecType const& vectorIn, std::vector<unsigned> const& indsIn) : vector(vectorIn), inds(indsIn){};

    double const& operator()(unsigned row) const{return vector(inds.at(row));};

    unsigned dimension(unsigned dim) const
    {
	if(dim>0)
	    return 1;
	else
	    return inds.size();
    }
    
private:
    VecType const& vector;
    std::vector<unsigned> const& inds;
};

template<typename MatType>
class MatrixBlock
{
public:
    MatrixBlock(MatType & matrixIn,
		unsigned startRowIn,
		unsigned startColIn,
		unsigned numRowsIn,
		unsigned numColsIn) : startRow(startRowIn), startCol(startColIn), numRows(numRowsIn), numCols(numColsIn), matrix(matrixIn)
    {
	assert(GetShape(matrix,0)>=startRow+numRows);
	assert(GetShape(matrix,1)>=startCol+numCols);
    };

    double& operator()(unsigned row, unsigned col) const{return matrix(startRow + row, startCol + col);};

    template<typename Derived>
    MatrixBlock& operator=(Eigen::DenseBase<Derived> const& otherMat)
    {
	assert(otherMat.rows()==numRows);
	assert(otherMat.cols()==numCols);
	
	for(int j=0; j<otherMat.cols(); ++j)
	{
	    for(int i=0; i<otherMat.rows(); ++i)
		matrix(i,j) = otherMat(i,j);
	}
	return *this;
    }
    
    unsigned rows() const{return numRows;};
    unsigned cols() const{return numCols;};
    
private:
    
    const unsigned startRow, startCol, numRows, numCols;
    MatType        & matrix;
};

/** 

Grab a particular column of a matrix and return as an instance of the "ColumnSlice" class. 

*/
template<typename MatType>
ColumnSlice<MatType> GetColumn(MatType const& matrix, unsigned col)
{
    return ColumnSlice<MatType>(matrix, col);
}

template<typename MatType>
VectorSlice<MatType> GetSlice(MatType const& matrix, std::vector<unsigned> const& inds)
{
    return VectorSlice<MatType>(matrix, inds);
}


template<typename MatType>
MatrixBlock<MatType> GetBlock(MatType & matrix, unsigned rowStart, unsigned colStart, unsigned numRows, unsigned numCols)
{
    return MatrixBlock<MatType>(matrix, rowStart, colStart, numRows, numCols);
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
	KernelEvaluator<StartDim, EndDim, KTypes...>::Evaluate(kernels, xs, ys, GetColumn(xs, xcol), GetColumn(ys,ycol), output);
    }

    template<typename MatrixType>
    static void Evaluate(std::tuple<KTypes...> const& kernels,
			 MatrixType            const& xs,
			 MatrixType            const& ys,
			 Eigen::VectorXd            & output)
    {
	output(StartDim) = std::get<StartDim>(kernels)( xs, ys );
        KernelEvaluator<StartDim+1, EndDim, KTypes...>::Evaluate(kernels, xs, ys, output);
	
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
    static void GetParams(std::tuple<KTypes...>     const& kernels,
			  Eigen::DenseCoeffsBase<Derived> & params,
	                  unsigned                         currInd =0)
    {
        auto kernel = std::get<StartDim>(kernels);
	unsigned numParams = kernel.GetNumParams();
	params.segment(currInd,numParams) = kernel.GetParams();
	
        KernelEvaluator<StartDim+1, EndDim, KTypes...>::GetParams(kernels, params, currInd+numParams);
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

    template<typename MatrixType>
    static double GetDeriv(std::tuple<KTypes...>     const& kernels,
			   MatrixType                const& x1,
			   MatrixType                const& x2,
			   unsigned                         wrt,
			   unsigned                         currParams, // the number of parameters in all previous kernels
	                   unsigned                       & kernelInd) // The current kernel of interest (same as the start dim template)
    {
        auto kernel = std::get<StartDim>(kernels);
	unsigned numParams = kernel.GetNumParams();
	if(wrt<currParams + numParams)
	{
	    Eigen::MatrixXd temp(1,1);
	    kernel.GetDerivative(x1, x2, wrt-currParams, temp);
	    kernelInd = StartDim;
	    return temp(0,0);
	}
	else
	{   
	    KernelEvaluator<StartDim+1, EndDim, KTypes...>::SetParams(kernels, x1,x2, wrt, currParams + numParams, kernelInd);
	}
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
    static void GetParams(std::tuple<KTypes...>     const& kernels,
			  Eigen::DenseCoeffsBase<Derived> & params,
	                  unsigned                         currInd = 0)
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
