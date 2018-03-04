#include "MUQ/Approximation/GaussianProcesses/CoregionalKernel.h"


CoregionalKernel::CoregionalKernel(unsigned                                        dim,
								 									 Eigen::MatrixXd                          const& Gamma,
								 								 	 std::vector<std::shared_ptr<KernelBase>> const& kernelsIn) : KernelBase(dim, Gamma.rows(), GetNumParams(kernelsIn)),
																																							kernels(kernelsIn)
{
	 // Make sure the matrix and kernels are the same size
	 assert(Gamma.cols()==kernelsIn.size());

	 // Compute the eigenvalue decomposition of the covariance Gamma to get the matrix A
	 Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver;
	 eigSolver.compute(Gamma, Eigen::ComputeEigenvectors);
	 A = eigSolver.eigenvectors()*eigSolver.eigenvalues().cwiseSqrt().asDiagonal();
};
