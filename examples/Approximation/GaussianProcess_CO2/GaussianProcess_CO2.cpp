#include <iostream>

#include "MUQ/Utilities/HDF5/H5Object.h"

#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

using namespace std;
using namespace muq::Utilities;
using namespace muq::Approximation;

int main()
{

    string dataFile = "data/MaunaLoaCO2.h5";

    //H5Object f = OpenFile(dataFile);

    //Eigen::VectorXd times          = f["/CO2/Times" ];
    //Eigen::VectorXd concentrations = f["/CO2/Concentrations" ];

    double k1_var = 10.0;
    double k1_length = 20;
    auto k1  = SquaredExpKernel(1, k1_var, k1_length);

    double k2_var = 1.0;
    double k2_period = 1.0;
    double k2_length = 2.0;
    auto k2 = PeriodicKernel(1, k2_var, k2_length, k2_period);

    double k3_var = 1.0;
    double k3_length = 20;
    double k3_nu = 3.0/2.0;
    auto k3 = MaternKernel(1, k3_var, k3_length, k3_nu);

    double k4_var = 1.0;
    double k4_length = 10;
    double k4_nu = 3.0/2.0;
    auto k4 = MaternKernel(1, k4_var, k4_length, k4_nu);

    // Combine the individual kernels
    auto k = k1 + k2*k3 + k4;


    ConstantMean mu(1, 1);
    GaussianProcess gp(mu, k);


    int numPts = 1000;
    Eigen::MatrixXd evalPts(1, numPts);
    evalPts.row(0) = Eigen::VectorXd::LinSpaced(numPts, 2000, 2020);

    Eigen::MatrixXd samp = gp.Sample(evalPts);

    std::cout << "Points,Concentration\n";
    for(int i=0; i<numPts; ++i)
        cout << evalPts(i) << "," << samp(i) << std::endl;

    return 0;
}
