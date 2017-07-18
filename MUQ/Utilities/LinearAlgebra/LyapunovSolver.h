#ifndef LYAPUNOVSOLVER_H
#define LYAPUNOVSOLVER_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <iostream>

namespace muq
{
namespace Utilities
{

    /** Solve the equation AX + XA^T + Q = 0
     */
    template<class ScalarType, int FixedRows=Eigen::Dynamic, int FixedCols=Eigen::Dynamic>
    class LyapunovSolver
    {
    public:

        typedef Eigen::Matrix<ScalarType, FixedRows, FixedCols> MatrixType;
        typedef Eigen::Matrix<std::complex<ScalarType>, FixedRows, FixedCols> ComplexMatrixType;
        
        void compute(MatrixType const& A, MatrixType const& C)
        {
            const int dim = A.rows();
            
            Eigen::ComplexSchur<MatrixType> schur;
            schur.compute(A);

            auto& Q = schur.matrixU();
            ComplexMatrixType S = schur.matrixT();

            ComplexMatrixType ctilde = Q.transpose() * C.template cast<std::complex< ScalarType >>() * Q;

            X.resize(dim,dim);
            ComputeFromSchur(ctilde,S,X);
        };

        MatrixType const& matrixX() const{return X;};

    private:

        void ComputeFromSchur(Eigen::Ref<const ComplexMatrixType> const& S,
                              Eigen::Ref<ComplexMatrixType>              ctilde,
                              Eigen::Ref<MatrixType>                     X)
        {
            const int size = X.rows();
            
            X(0,0) = -(ctilde(0,0)/(S(0,0)+ std::conj(S(0,0)))).real();

            if(size==1)
                return;
            
            ComplexMatrixType tempX = -ctilde.block(1,0,size-1,1) - X(0,0)*S.block(0,1,1,size-1).transpose();
            (S.block(1,1,size-1,size-1).transpose() + S(0,0)*ComplexMatrixType::Identity(size-1,size-1)).template triangularView<Eigen::Lower>().solveInPlace(tempX);

            X.block(1,0,size-1,1) = tempX.real();
            X.block(0,1,1,size-1) = X.block(1,0,size-1,1).transpose();

            // Recursively call this function on the next block
            ctilde.block(1,1,size-1,size-1) += S.block(0,1,1,size-1).transpose()*X.block(0,1,1,size-1) + X.block(1,0,size-1,1)*S.block(0,1,1,size-1);
                    
            ComputeFromSchur(S.block(1,1,size-1,size-1), ctilde , X.block(1,1,size-1,size-1));
            
        }
    
        MatrixType X;
        
    };
}
}


#endif
