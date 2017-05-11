#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

using namespace muq::Approximation;


double muq::Approximation::nlopt_obj(unsigned n, const double *x, double *nlopt_grad, void *opt_info)
{
    Eigen::VectorXd grad;
    Eigen::Map<const Eigen::VectorXd> params(x,n);
    double logLikely;
    
    OptInfo* info = (OptInfo*) opt_info;
    
    info->gp->covKernel->SetParams(params);
    
    logLikely = info->gp->EvaluateMarginalLikelihood(*info->xs, *info->vals, grad, nlopt_grad);	    
    
    if (nlopt_grad)
    {
	for(int i=0; i<n; ++i)
	    nlopt_grad[i] = grad(i);
    }
    
    return logLikely;
}


GaussianProcess::GaussianProcess(std::shared_ptr<MeanFunctionBase> meanIn,
				 std::shared_ptr<KernelBase>       covKernelIn) :
    mean(meanIn),
    covKernel(covKernelIn),
    inputDim(covKernel->inputDim),
    coDim(covKernel->coDim)
{};

void GaussianProcess::Fit(Eigen::MatrixXd const& xs,
			  Eigen::MatrixXd const& vals)
{
    assert(xs.cols()==vals.cols());
    
    // Store the training data
    trainLocs = xs;
    trainVals = vals;
    
    // Further setup...
    Fit();
}

void GaussianProcess::Fit()
{
    // Evaluate the priorMean
    trainMean = mean->Evaluate(trainLocs);
    
    // Store the covariance information
    Eigen::MatrixXd cov(trainLocs.cols(), trainLocs.cols());
    covKernel->FillCovariance(trainLocs, cov);
    covSolver = cov.llt();
}


void GaussianProcess::Optimize()
{
    
    OptInfo info;
    info.xs = &trainLocs;
    info.vals = &trainVals;
    info.gp = this;
    
    Eigen::VectorXd params = covKernel->GetParams();

    const Eigen::MatrixXd bounds = covKernel->GetParamBounds();
    Eigen::VectorXd lbs = bounds.row(0);
    Eigen::VectorXd ubs = bounds.row(1);
	    

    nlopt_opt opt;

    opt = nlopt_create(NLOPT_LD_LBFGS, params.rows()); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lbs.data());
    nlopt_set_upper_bounds(opt, ubs.data());
    nlopt_set_vector_storage(opt, 100);
    nlopt_set_max_objective(opt, muq::Approximation::nlopt_obj, (void*) &info);
	    	    
    double maxLikely;
    nlopt_optimize(opt, params.data(), &maxLikely);

    covKernel->SetParams(params);
	    
    nlopt_destroy(opt);

    // Now that we've optimized the hyperparameters, get set up for predictions
    Fit(trainLocs, trainVals);

}


GaussianInformation GaussianProcess::Predict(Eigen::MatrixXd const& newLocs)
{

    // Construct the cross covariance
    Eigen::MatrixXd crossCov(newLocs.cols(), trainLocs.cols());
    covKernel->FillCovariance(newLocs, trainLocs, crossCov);

    GaussianInformation output;
		    
    // Solve for the posterior mean
    Eigen::Map<Eigen::VectorXd> valMap(trainVals.data(), trainVals.rows()*trainVals.cols());
    Eigen::Map<Eigen::VectorXd> priorMeanMap(trainMean.data(), trainMean.rows()*trainMean.cols());

    output.mean.resize(trainVals.rows(), newLocs.cols());
    Eigen::Map<Eigen::VectorXd> postMean(output.mean.data(), trainVals.rows()*newLocs.cols());
	    
    postMean = crossCov * covSolver.solve(valMap - priorMeanMap);
    output.mean += mean->Evaluate(newLocs);
	    
    // Compute the prior covariance
    Eigen::MatrixXd priorCov(newLocs.cols(), newLocs.cols());
    covKernel->FillCovariance(newLocs,priorCov);
	    
    // Solve for the posterior covariance
    output.covariance = priorCov - crossCov * covSolver.solve(crossCov.transpose());

    return output;
};


// Evaluates the log marginal likelihood needed when fitting hyperparameters
double  GaussianProcess::EvaluateMarginalLikelihood(Eigen::MatrixXd const& xs,
						    Eigen::MatrixXd const& vals,
						    Eigen::VectorXd      & grad,
						    bool                   computeGrad)
{

    int numXs = xs.cols();

    // build the covariance matrix between the input points
    Eigen::MatrixXd cov(coDim * numXs, coDim * numXs);
    covKernel->FillCovariance(xs, cov);

    // Compute the cholesky decomposition of the covariance
    Eigen::LLT<Eigen::MatrixXd> chol = cov.llt();

    // Compute the log determinant of the covariance
    const Eigen::TriangularView<const Eigen::MatrixXd, Eigen::Lower> L = chol.matrixL();
    double logDet = 0.0;
    for(int i=0; i<L.rows(); ++i)
	logDet += 2.0*log(L(i,i));
	    
    // Make the mean prediction and compute the difference with observations
    Eigen::MatrixXd obsDiff = vals - mean->Evaluate(xs);
    Eigen::Map<Eigen::VectorXd> unwrappedDiff(obsDiff.data(), coDim*numXs);

    // Compute the log likelihood
    Eigen::VectorXd alpha = chol.solve(unwrappedDiff);
    double logLikely = -0.5 * unwrappedDiff.transpose() * alpha - 0.5*logDet - 0.5*numXs*log(2.0*pi);

    if(computeGrad)
    {
	// Compute the gradient of the log likelihood
	const int numParams = covKernel->numParams;
	grad.resize(numParams);
	Eigen::MatrixXd derivMat(coDim * numXs, coDim * numXs);
	for(int i=0; i<numParams; ++i)
	{
	    covKernel->FillDerivativeMatrix(xs, i, derivMat);
	    grad(i) = 0.5*(alpha*alpha.transpose()*derivMat - chol.solve(derivMat)).trace();
	}
    }

    return logLikely;	    
};
