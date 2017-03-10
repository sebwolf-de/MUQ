#ifndef GAUSSIANPROCESS_H_
#define GAUSSIANPROCESS_H_

#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"

#include <Eigen/Core>
#include <set>

#include <nlopt.h>


namespace muq
{

    namespace Approximation
    {
	

    template<typename MeanType, typename KernelType>
    class GaussianProcess;

    template<typename MeanType, typename KernelType>
    struct OptInfo
    {
	const Eigen::MatrixXd *xs;
	const Eigen::MatrixXd *vals;

	GaussianProcess<MeanType, KernelType> *gp;
    };
    
    template<typename MeanType, typename KernelType>
    double nlopt_obj(unsigned n, const double *x, double *nlopt_grad, void *opt_info)
    {
	Eigen::VectorXd grad;
	Eigen::Map<const Eigen::VectorXd> params(x,n);
	double logLikely;

	OptInfo<MeanType,KernelType>* info = (OptInfo<MeanType,KernelType>*) opt_info;

	info->gp->covKernel.SetParams(params);
	
	logLikely = info->gp->EvaluateLikelihood(*info->xs, *info->vals, grad, nlopt_grad);	    

	if (nlopt_grad)
	{
	    for(int i=0; i<n; ++i)
		nlopt_grad[i] = grad(i);
	}

	return logLikely;
    }

	
    
    class ConstantMean
    {

    public:
        ConstantMean(int coDimIn) : coDim(coDimIn){};

	
	Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs)
	{
	    return Eigen::MatrixXd::Zero(coDim, xs.cols());
	}

    private:

	const int coDim;
	    
    };


    class GaussianInformation
    {
    public:
	Eigen::MatrixXd mean;
	Eigen::MatrixXd covariance;
    };
    
    template<typename MeanType, typename KernelType>
    class GaussianProcess
    {

    public:
	GaussianProcess(MeanType   const& priorMeanIn,
			KernelType const& covKernelIn) :
	    priorMean(priorMeanIn),
	    covKernel(covKernelIn),
	    inputDim(covKernel.GetDim()),
	    coDim(covKernel.GetCodim())
        {};


	    
	void Fit(Eigen::MatrixXd const& xs,
		 Eigen::MatrixXd const& vals)
	{
	    assert(xs.cols()==vals.cols());

	    // Store the training data
	    trainLocs = xs;
	    trainVals = vals;

	    // Further setup...
	    Fit();
	}

	void Fit()
	{
	    // Evaluate the priorMean
	    trainMean = priorMean.Evaluate(trainLocs);

	    // Store the covariance information
	    Eigen::MatrixXd cov(trainLocs.cols(), trainLocs.cols());
	    covKernel.BuildCovariance(trainLocs, cov);
	    covSolver = cov.llt();
	}

	void Optimize()
	{

	    OptInfo<MeanType, KernelType> info;
	    info.xs = &trainLocs;
	    info.vals = &trainVals;
	    info.gp = this;

	    Eigen::VectorXd params = covKernel.GetParams();

	    const Eigen::MatrixXd bounds = covKernel.GetParamBounds();
	    Eigen::VectorXd lbs = bounds.row(0);
	    Eigen::VectorXd ubs = bounds.row(1);
	    

	    nlopt_opt opt;

	    opt = nlopt_create(NLOPT_LD_LBFGS, params.rows()); /* algorithm and dimensionality */
	    nlopt_set_lower_bounds(opt, lbs.data());
	    nlopt_set_upper_bounds(opt, ubs.data());
	    nlopt_set_vector_storage(opt, 100);
	    nlopt_set_max_objective(opt, nlopt_obj<MeanType,KernelType>, (void*) &info);
	    	    
	    double maxLikely;
	    nlopt_optimize(opt, params.data(), &maxLikely);

	    covKernel.SetParams(params);
	    
	    nlopt_destroy(opt);

	    // Now that we've optimized the hyperparameters, get set up for predictions
	    Fit(trainLocs, trainVals);

	}

	GaussianInformation Predict(Eigen::MatrixXd const& newLocs)
	{

	    // Construct the cross covariance
	    Eigen::MatrixXd crossCov(newLocs.cols(), trainLocs.cols());
	    covKernel.BuildCovariance(newLocs, trainLocs, crossCov);

	    GaussianInformation output;

	    std::cout << "Here 0" << std::endl;
		    
	    // Solve for the posterior mean
	    Eigen::Map<Eigen::VectorXd> valMap(trainVals.data(), trainVals.rows()*trainVals.cols());
	    Eigen::Map<Eigen::VectorXd> priorMeanMap(trainMean.data(), trainMean.rows()*trainMean.cols());

	    output.mean.resize(trainVals.rows(), newLocs.cols());
	    Eigen::Map<Eigen::VectorXd> postMean(output.mean.data(), trainVals.rows()*newLocs.cols());

	    std::cout << "Here 1" << std::endl;
	    
	    postMean = crossCov * covSolver.solve(valMap - priorMeanMap);
            output.mean += priorMean.Evaluate(newLocs);
	    
	    std::cout << "Here 2" << std::endl;
	    
	    // Compute the prior covariance
	    Eigen::MatrixXd priorCov(newLocs.cols(), newLocs.cols());
	    covKernel.BuildCovariance(newLocs,priorCov);
	    
	    // Solve for the posterior covariance
	    output.covariance = priorCov - crossCov * covSolver.solve(crossCov.transpose());

	    return output;
	};
	
	// Evaluates the log marginal likelihood needed when fitting hyperparameters
	double EvaluateLikelihood(Eigen::MatrixXd const& xs,
				  Eigen::MatrixXd const& vals,
				  Eigen::VectorXd      & grad,
	                          bool                   computeGrad = true)
	{

	    int numXs = xs.cols();

	    // build the covariance matrix between the input points
	    Eigen::MatrixXd cov(coDim * numXs, coDim * numXs);
	    covKernel.BuildCovariance(xs, cov);

	    // Compute the cholesky decomposition of the covariance
	    Eigen::LLT<Eigen::MatrixXd> chol = cov.llt();

	    // Compute the log determinant of the covariance
	    auto L = chol.matrixL();
	    double logDet = 0.0;
	    for(int i=0; i<L.rows(); ++i)
		logDet += 2.0*log(L(i,i));
	    
	    // Make the mean prediction and compute the difference with observations
	    Eigen::MatrixXd obsDiff = vals - priorMean.Evaluate(xs);
	    Eigen::Map<Eigen::VectorXd> unwrappedDiff(obsDiff.data(), coDim*numXs);

	    // Compute the log likelihood
	    Eigen::VectorXd alpha = chol.solve(unwrappedDiff);
	    double logLikely = -0.5 * unwrappedDiff.transpose() * alpha - 0.5*logDet - 0.5*numXs*log(2.0*pi);

	    if(computeGrad)
	    {
		// Compute the gradient of the log likelihood
		const int numParams = covKernel.GetNumParams();
		grad.resize(numParams);
		Eigen::MatrixXd derivMat(coDim * numXs, coDim * numXs);
		for(int i=0; i<numParams; ++i)
		{
		    covKernel.BuildDerivativeMatrix(xs, i, derivMat);
		    grad(i) = 0.5*(alpha*alpha.transpose()*derivMat - chol.solve(derivMat)).trace();
		}
	    }

	    return logLikely;	    
	};

	
	Eigen::MatrixXd EvaluateMean(Eigen::MatrixXd const& xs);


	MeanType   priorMean;
	KernelType covKernel;

    private:

	Eigen::MatrixXd trainMean;
	Eigen::MatrixXd trainLocs;
	Eigen::MatrixXd trainVals;
	
	Eigen::LLT<Eigen::MatrixXd> covSolver;
	
	const int inputDim;
	const int coDim;

	const double pi = 4.0 * atan(1.0); //boost::math::constants::pi<double>();
	
    };// class GaussianProcess


    template<typename MeanType, typename KernelType>
	GaussianProcess<MeanType,KernelType> ConstructGP(MeanType const& mean,
							 KernelType const& kernel)
    {
	return GaussianProcess<MeanType,KernelType>(mean,kernel);
    }

    } // namespace Approximation
} // namespace muq

#endif // #ifndef GAUSSIANPROCESS_H_
