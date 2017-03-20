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
	
    class GaussianProcess;

    
    struct OptInfo
    {
	const Eigen::MatrixXd *xs;
	const Eigen::MatrixXd *vals;

	GaussianProcess *gp;
    };
    
    double nlopt_obj(unsigned n, const double *x, double *nlopt_grad, void *opt_info);

    class MeanFunctionBase
    {
    public:
        MeanFunctionBase(unsigned dimIn,
	                 unsigned coDimIn) : dim(dimIn), coDim(coDimIn){}
	
	virtual Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs) const = 0;

	virtual std::shared_ptr<MeanFunctionBase> Clone() const = 0;

    protected:
	unsigned dim, coDim;
    };
    
    class ConstantMean : public MeanFunctionBase
    {

    public:
        ConstantMean(unsigned dim, unsigned coDim) : MeanFunctionBase(dim,coDim){};

        virtual std::shared_ptr<MeanFunctionBase> Clone() const override
	{
	    return std::make_shared<ConstantMean>(*this);
	}
		
	virtual Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs) const override
	{
	    return Eigen::MatrixXd::Zero(coDim, xs.cols());
	}	    
    };


    class GaussianInformation
    {
    public:
	Eigen::MatrixXd mean;
	Eigen::MatrixXd covariance;
    };
    
    class GaussianProcess
    {

    public:
    GaussianProcess(std::shared_ptr<MeanFunctionBase> meanIn,
		    std::shared_ptr<KernelBase>       covKernelIn);


	    
	void Fit(Eigen::MatrixXd const& xs,
		 Eigen::MatrixXd const& vals);

	void Fit();

	void Optimize();

	GaussianInformation Predict(Eigen::MatrixXd const& newLocs);
	
	// Evaluates the log marginal likelihood needed when fitting hyperparameters
	double EvaluateMarginalLikelihood(Eigen::MatrixXd const& xs,
					  Eigen::MatrixXd const& vals,
					  Eigen::VectorXd      & grad,
					  bool                   computeGrad = true);

	
	Eigen::MatrixXd EvaluateMean(Eigen::MatrixXd const& xs);

	std::shared_ptr<MeanFunctionBase> mean;
	std::shared_ptr<KernelBase>       covKernel;

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
    GaussianProcess ConstructGP(MeanType const& mean,
				KernelType const& kernel)
    {
	return GaussianProcess(mean.Clone(),kernel.Clone());
    }

    
    } // namespace Approximation
} // namespace muq

#endif // #ifndef GAUSSIANPROCESS_H_
