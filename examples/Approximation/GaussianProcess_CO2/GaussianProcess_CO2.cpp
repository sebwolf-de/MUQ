#include <iostream>

#include "MUQ/Utilities/HDF5/H5Object.h"

#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

using namespace std;
using namespace muq::Utilities;
using namespace muq::Approximation;

int main()
{
    // Open and read in data
    string dataFile = "data/MaunaLoaCO2.h5";
    H5Object f = OpenFile(dataFile);

    Eigen::VectorXd times          = f["/Weekly/Dates" ];
    Eigen::VectorXd concentrations = f["/Weekly/Concentrations" ];

    double k1_var = 360;
    double k1_length = 70;
    auto k1  = SquaredExpKernel(1, k1_var, k1_length);

    double k2_var = 1.0;
    double k2_period = 1.0;
    double k2_length = 2.0;
    auto k2 = PeriodicKernel(1, k2_var, k2_length, k2_period);

    double k3_var = 5.0;
    double k3_length = 90;
    double k3_nu = 3.0/2.0;
    auto k3 = MaternKernel(1, k3_var, k3_length, k3_nu);

    double k4_var = 0.5;
    double k4_length = 2.0;
    double k4_nu = 3.0/2.0;
    auto k4 = MaternKernel(1, k4_var, k4_length, k4_nu);

    // Combine the individual kernels
    auto k = k1 + k2*k3 + k4;

    ConstantMean mu(1, 1);
    GaussianProcess gp(mu, k);


    // Define pediction points
    int numPts = 500;
    Eigen::MatrixXd evalPts(1, numPts);
    evalPts.row(0) = Eigen::VectorXd::LinSpaced(numPts, 2000, 2050);

    gp.Condition(times.transpose(), concentrations.transpose(), 1e-2);
    Eigen::MatrixXd postMean, postVar;
    std::tie(postMean,postVar) = gp.Predict(evalPts, GaussianProcess::DiagonalCov);

    // Write results to new file
    string writeFile = "results/CO2_Prediction.h5";
    H5Object fout = OpenFile(writeFile);

    fout["/Predict/Dates"] = evalPts;
    fout["/Predict/Concentrations"] = postMean;
    fout["/Predict/ConcentrationVariance"] = postVar;
    return 0;
}
