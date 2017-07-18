#include "MUQ/Utilities/LinearAlgebra/LyapunovSolver.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

using namespace muq::Utilities;


TEST(Utilities_LyapunovSolver, Diagonal)
{
    const int dim = 10;
    Eigen::MatrixXd A = -Eigen::MatrixXd::Identity(dim,dim);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(dim,dim);


    LyapunovSolver<double> solver;
    solver.compute(A,Q);

    Eigen::MatrixXd const& X = solver.matrixX();

    for(int i=0; i<dim; ++i)
    {
        EXPECT_DOUBLE_EQ(0.5, X(i,i) );

        for(int j=i+1; j<dim; ++j)
        {
            EXPECT_DOUBLE_EQ(0.0, X(i,j));
            EXPECT_DOUBLE_EQ(0.0, X(j,i));
        }
    }

}
