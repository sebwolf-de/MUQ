#include "MUQ/Utilities/LinearAlgebra/ConcatenateOperator.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Utilities;

ConcatenateOperator::ConcatenateOperator(std::shared_ptr<LinearOperator> Ain,
                                         std::shared_ptr<LinearOperator> Bin,
                                         const int                       rowColIn) : LinearOperator(GetRows(Ain,Bin, rowColIn), GetCols(Ain,Bin,rowColIn)), A(Ain), B(Bin), rowCol(rowColIn)
{
    CheckSizes();
};


/** Apply the linear operator to a vector */
Eigen::MatrixXd ConcatenateOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    Eigen::MatrixXd output(nrows, x.cols());
    
    if(rowCol==0){
        output.topRows(A->rows()) = A->Apply(x);
        output.bottomRows(B->rows()) = B->Apply(x);
    }else{
        output = A->Apply(x.topRows(A->cols()));
        output += B->Apply(x.bottomRows(B->cols()));
    }
    
    return output;
}
  

/** Apply the transpose of the linear operator to a vector. */
Eigen::MatrixXd ConcatenateOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    Eigen::MatrixXd output(ncols, x.cols());
    
    if(rowCol==0){
        output = A->ApplyTranspose( x.topRows(A->rows()) );
        output += B->ApplyTranspose( x.bottomRows(B->rows()) );
    }else{
        output.topRows(A->cols()) = A->ApplyTranspose( x.topRows(A->rows()) );
        output.bottomRows(B->cols()) = B->ApplyTranspose(x.bottomRows(B->rows()));
    }
    
    return output;

}
    

 std::shared_ptr<ConcatenateOperator> ConcatenateOperator::VStack(std::shared_ptr<LinearOperator> Ain,
                                                                  std::shared_ptr<LinearOperator> Bin)
 {
     return std::make_shared<ConcatenateOperator>(Ain,Bin,0);
 }

std::shared_ptr<ConcatenateOperator> ConcatenateOperator::HStack(std::shared_ptr<LinearOperator> Ain,
                                                                 std::shared_ptr<LinearOperator> Bin)
{
    return std::make_shared<ConcatenateOperator>(Ain,Bin,1);
}


Eigen::MatrixXd ConcatenateOperator::GetMatrix()
{

    Eigen::MatrixXd output(nrows, ncols);
    if(rowCol==0){
        output.topRows(A->rows()) = A->GetMatrix();
        output.bottomRows(A->rows()) = B->GetMatrix();
    }else{
        output.leftCols(A->cols()) = A->GetMatrix();
        output.rightCols(B->cols()) = B->GetMatrix();
    }
    return output;
}

int ConcatenateOperator::GetRows(std::shared_ptr<LinearOperator> Ain,
                                 std::shared_ptr<LinearOperator> Bin,
                                 const int                       rowColIn)
{
    if(rowColIn==0){
        return Ain->rows() + Bin->rows();
    }else{
        return Ain->rows();
    }
}


int ConcatenateOperator::GetCols(std::shared_ptr<LinearOperator> Ain,
                                 std::shared_ptr<LinearOperator> Bin,
                                 const int                       rowColIn)
{
    if(rowColIn==0){
        return Ain->cols();
    }else{
        return Ain->cols() + Bin->cols();
    }
}


void ConcatenateOperator::CheckSizes()
{

    if(rowCol==0){
        if(A->cols() != B->cols())
            throw muq::WrongSizeError("In ConcatenateOperator: Cannot vertically stack operators with different number of columns.  Matrix A has " + std::to_string(A->cols()) + " columns but matrix B has " + std::to_string(B->cols()) + " columns.");
        
    }else{
        if(A->rows() != B->rows())
            throw muq::WrongSizeError("In ConcatenateOperator: Cannot horizontally stack operators with different number of rows.  Matrix A has " + std::to_string(A->rows()) + " rows but matrix B has " + std::to_string(B->rows()) + " rows.");

    }
}
