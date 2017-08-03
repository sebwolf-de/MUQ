#ifndef BLOCKDATASET_H
#define BLOCKDATASET_H

#include "MUQ/Utilities/HDF5/HDF5File.h"

namespace muq
{
namespace Utilities
{
    
    class BlockDataset
    {
	
    public:

	BlockDataset(std::string               const& path_,
		     std::shared_ptr<HDF5File>        file_,
	             const int                        startRow_,
		     const int                        startCol_,
		     const int                        numRows_,
		     const int                        numCols_) : path(path_),
	                                                          file(file_),
	                                                          startRow(startRow_),
	                                                          startCol(startCol_),
	                                                          numRows(numRows_),
	                                                          numCols(numCols_){};
	
	template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
	BlockDataset& operator=(ScalarType val)
	{
	    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> temp = val*Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Ones(numRows, numCols);

	    file->WritePartialMatrix(path, temp, startRow, startCol);
	    
	    return *this;
	}

	template<typename MatrixType, typename = typename std::enable_if< !std::is_arithmetic<MatrixType>::value, MatrixType>::type>
	BlockDataset& operator=(MatrixType const& val)
	{
	    assert(val.rows()==numRows);
	    assert(val.cols()==numCols);
	    
	    file->WritePartialMatrix(path, val, startRow, startCol); 
	    return *this;
	}

	template<typename ScalarType, int rows, int cols>
	operator Eigen::Matrix<ScalarType,rows,cols>()
	{
	    return eval<ScalarType,rows,cols>();
	}

	template<typename ScalarType=double, int rows=Eigen::Dynamic, int cols=Eigen::Dynamic>
	Eigen::Matrix<ScalarType,rows,cols> eval()
	{
	    return file->ReadPartialMatrix<ScalarType, rows, cols>(path, startRow, startCol, numRows, numCols);
	}
	
	template<typename ScalarType, int rows, int cols>
	operator ScalarType()
	{
	    assert(numRows==1);
	    assert(numCols==1);
	    
	    return file->ReadPartialMatrix<ScalarType, rows, cols>(path, startRow, startCol, numRows, numCols)(0);
	}

	
	
    private:

	// Path to the original dataset
	const std::string path;

	// The HDF5File
	std::shared_ptr<HDF5File> file;

	// Integers defining the slice of the dataset
	const int startRow;
	const int startCol;
	const int numRows;
	const int numCols;
	
    }; // class BlockDataset

    
} // namespace Utilities
} // namespace muq


#endif
