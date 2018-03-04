#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

using namespace muq::Approximation;

Eigen::MatrixXd KernelBase::Evaluate(Eigen::VectorXd const& x1,
																     Eigen::VectorXd const& x2) const
{
	Eigen::MatrixXd output(coDim, coDim);
	FillBlock(x1(dimInds), x2(dimInds), cachedParams, output);
	return output;
}

Eigen::MatrixXd KernelBase::BuildCovariance(Eigen::MatrixXd const& x) const
{
	Eigen::MatrixXd output(coDim*x.cols(), coDim*x.cols());
	FillCovariance(x,output);
	return output;
}

Eigen::MatrixXd KernelBase::BuildCovariance(Eigen::MatrixXd const& x1,
	                                          Eigen::MatrixXd const& x2) const
{
	Eigen::MatrixXd output(coDim*x1.cols(), coDim*x2.cols());
	FillCovariance(x1,x2,output);
	return output;
}

void KernelBase::FillCovariance(Eigen::MatrixXd      const& x,
                                Eigen::Ref<Eigen::MatrixXd> output) const
{
	assert(output.rows()==x.cols()*coDim);
	assert(output.cols()==x.cols()*coDim);

	const int numPts = x.cols();

	#pragma omp parallel for
	for(int i=0; i<numPts; ++i){

		for(int j=0; j<i; ++j){
			FillBlock(x.col(i)(dimInds), x.col(j)(dimInds), cachedParams, output.block(i*coDim, j*coDim, coDim, coDim));
			output.block(j*coDim, i*coDim, coDim, coDim) = output.block(i*coDim, j*coDim, coDim, coDim).transpose();
		}

		FillBlock(x.col(i)(dimInds), x.col(i)(dimInds), cachedParams, output.block(i*coDim,i*coDim,coDim,coDim));
	}

}

void KernelBase::FillCovariance(Eigen::MatrixXd      const& x1,
																Eigen::MatrixXd      const& x2,
															  Eigen::Ref<Eigen::MatrixXd> output) const
{
	int numRows = x1.cols();
	int numCols = x2.cols();

	assert(output.rows()==numRows*coDim);
	assert(output.cols()==numCols*coDim);

	#pragma omp parallel for
	for(int i=0; i<numRows; ++i){
		for(int j=0; j<numCols; ++j)
			FillBlock(x1.col(i)(dimInds), x2.col(j)(dimInds), cachedParams, output.block(i*coDim, j*coDim, coDim, coDim));
	}
}
